#include "./utils.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/getpot.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

static void input (GetPot & , EquationSystems & );
static void initial_structure (EquationSystems & , const std::string & );
static void initial_pihna (EquationSystems & , const std::string & );
static void assemble_pihna (EquationSystems & , const std::string & );

extern PerfLog plog;

void pihna (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);
  GetPot in("input.dat");

  input(in, es);

  TransientLinearImplicitSystem & model =
    es.add_system<TransientLinearImplicitSystem>("PIHNA");
  model.add_variable("n", FIRST, LAGRANGE);
  model.add_variable("c", FIRST, LAGRANGE);
  model.add_variable("h", FIRST, LAGRANGE);
  model.add_variable("v", FIRST, LAGRANGE);
  model.add_variable("a", FIRST, LAGRANGE);
  model.attach_assemble_function(assemble_pihna);
  model.attach_init_function(initial_pihna);

  ExplicitSystem & ustruct =
    es.add_system<ExplicitSystem>("uStructure");
  ustruct.add_variable("HU", CONSTANT, MONOMIAL);
  ustruct.add_variable("RT", CONSTANT, MONOMIAL);
  ustruct.attach_init_function(initial_structure);

  GmshIO(mesh).read(es.parameters.get<std::string>("input_GMSH"));
  mesh.prepare_for_use();
  mesh.print_info();
  GmshIO(mesh).write(es.parameters.get<std::string>("output_GMSH"));
  es.init();
  es.print_info();

  const std::string ex2_filename =
    es.parameters.get<std::string>("output_EXODUS");

  ExodusII_IO ex2(mesh);
  ex2.write_equation_systems(ex2_filename, es);
  ex2.append(true);

  const int output_step = es.parameters.get<int>("output_step");
  const int n_t_step = es.parameters.get<int>("time_step_number");
  for (int t=1; t<=n_t_step; t++)
    {
      es.parameters.set<Real>("time") +=
      es.parameters.get<Real>("time_step");
      // increment this time counter
      model.time = es.parameters.get<Real>("time");
      const Real current_time = model.time;

      libMesh::out << " Solving time increment: " << t
                   << " (time=" << current_time <<  ") ..." << std::endl;

      // copy the previously-current solution into the old solution
      *(model.old_local_solution) = *(model.current_local_solution);
      // now solve the AD progression model
      model.solve();

      if (0 == t%output_step)
        ex2.write_timestep(ex2_filename, es, t, current_time);
    }

  // ...done
}

void input (GetPot & in, EquationSystems & es)
{
  std::string name;

  // ...done
}

void initial_structure (EquationSystems & es,
                        const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "uStructure");

  const MeshBase& mesh = es.get_mesh();

  ExplicitSystem & system =
    es.get_system<ExplicitSystem>("uStructure");
  libmesh_assert_equal_to(system.n_vars(), 2);

  std::ifstream fin(es.parameters.get<std::string>("input_elemental"));

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      Real HU_, RT_;
      fin >> HU_ >> RT_;

      std::vector<std::vector<dof_id_type>> dof_indices_T_var(2);

      system.get_dof_map().dof_indices(elem, dof_indices_T_var[0], 0);
      system.solution->set(dof_indices_T_var[0][0], HU_);
      system.get_dof_map().dof_indices(elem, dof_indices_T_var[1], 1);
      system.solution->set(dof_indices_T_var[1][0], RT_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void initial_pihna (EquationSystems & es,
                    const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "PIHNA");

  const MeshBase& mesh = es.get_mesh();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PIHNA");
  libmesh_assert_equal_to(system.n_vars(), 5);

  es.parameters.set<Real> ("time") =
  system.time = 0.0;

  std::ifstream fin(es.parameters.get<std::string>("input_nodal"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real n_, c_, h_, v_, a_;
      fin >> n_ >> c_ >> h_ >> v_ >> a_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) ,
                                   node->dof_number(system.number(), 3, 0) ,
                                   node->dof_number(system.number(), 4, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );
      libmesh_assert( node->n_comp(system.number(), 3) == 1 );
      libmesh_assert( node->n_comp(system.number(), 4) == 1 );

      system.solution->set(idof[0], n_);
      system.solution->set(idof[1], c_);
      system.solution->set(idof[2], h_);
      system.solution->set(idof[3], v_);
      system.solution->set(idof[4], a_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void assemble_pihna (EquationSystems & es,
                     const std::string & system_name)
{
  libmesh_ignore(es, system_name);
  libmesh_assert_equal_to(system_name, "PIHNA");
  libmesh_assert_equal_to(system.n_vars(), 5);

  // ...done
}
