#include "./utils.h"

static void input (const std::string & , EquationSystems & );
static void initial_structure (EquationSystems & , const std::string & );
static void initial_pihna (EquationSystems & , const std::string & );
static void assemble_pihna (EquationSystems & , const std::string & );
static void check_solution (EquationSystems & );

extern PerfLog plog;

void pihna (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);

  input("input.dat", es);

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
      // update the simulation time
      es.parameters.set<Real>("time") += es.parameters.get<Real>("time_step");
      model.time = es.parameters.get<Real>("time");

      libMesh::out << " Solving time increment: " << t
                   << " (time=" << model.time <<  ") ..." << std::endl;

      // copy the previously-current solution into the old solution
      *(model.old_local_solution) = *(model.current_local_solution);
      // now solve the AD progression model
      model.solve();

      check_solution(es);

      if (0 == t%output_step)
        ex2.write_timestep(ex2_filename, es, t, model.time);
    }

  // ...done
}

void input (const std::string & file_name, EquationSystems & es)
{
  GetPot in(file_name);

  std::string name;

  const std::string DIR_default = date_time_to_string(date_now(), "%Y%m%d_%H%M%S");
  name = "directory";
  const std::string DIR = in(name, DIR_default) + "/";

  // create a directory to store simulation results
  if (0==global_processor_id())
    std::system(std::string("mkdir "+DIR).c_str());
  // create a copy of the input file containing all model parameters
  if (0==global_processor_id())
    std::system(std::string("cp "+file_name+" "+DIR+file_name).c_str());

  name = "input_GMSH";
  es.parameters.set<std::string>(name) = in(name, "input.msh");
  //
  name = "output_GMSH";
  es.parameters.set<std::string>(name) = DIR + in(name, "output.msh");
  //
  name = "input_nodal";
  es.parameters.set<std::string>(name) = in(name, "input.nodal");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());
  //
  name = "input_elemental";
  es.parameters.set<std::string>(name) = in(name, "input.elemental");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());
  //
  name = "output_EXODUS";
  es.parameters.set<std::string>(name) = DIR + in(name, "output.ex2");

  es.parameters.set<Real>("time") = 0.0;

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = in(name, 1);
  name = "output_step";
  es.parameters.set<int>(name) = in(name, 1);

  name = "cells_min_capacity";
  es.parameters.set<Real>(name) = in(name, 0.0);
  name = "cells_max_capacity";
  es.parameters.set<Real>(name) = in(name, 1.0);
  name = "cells_max_capacity/exponent";
  es.parameters.set<Real>(name) = in(name, 1.0);
  name = "cytokines_max_capacity";
  es.parameters.set<Real>(name) = in(name, 1.0);

  // parameters for the species: n
  {
    name = "necrosis/c";       es.parameters.set<Real>(name) = in(name, 0.);
    name = "necrosis/h";       es.parameters.set<Real>(name) = in(name, 0.);
    name = "necrosis/v";       es.parameters.set<Real>(name) = in(name, 0.);
  }

  // parameters for the species: c & h
  {
    name = "diffuse/c";        es.parameters.set<Real>(name) = in(name, 0.);
    name = "taxis/c";          es.parameters.set<Real>(name) = in(name, 0.);
    name = "diffuse/h";        es.parameters.set<Real>(name) = in(name, 0.);
    name = "taxis/h";          es.parameters.set<Real>(name) = in(name, 0.);
    name = "produce/c";        es.parameters.set<Real>(name) = in(name, 0.);
    name = "switch/c/to/h";    es.parameters.set<Real>(name) = in(name, 0.);
    name = "switch/h/to/c";    es.parameters.set<Real>(name) = in(name, 0.);
    name = "switch/h/to/n";    es.parameters.set<Real>(name) = in(name, 0.);
  }

  // parameters for the species: v
  {
    name = "diffuse/v";        es.parameters.set<Real>(name) = in(name, 0.);
    name = "taxis/v";          es.parameters.set<Real>(name) = in(name, 0.);
    name = "produce/v";        es.parameters.set<Real>(name) = in(name, 0.);
  }

  // parameters for the species: a
  {
    name = "secrete/a/from/c"; es.parameters.set<Real>(name) = in(name, 0.);
    name = "secrete/a/from/h"; es.parameters.set<Real>(name) = in(name, 0.);
    name = "uptake/a/from/v";  es.parameters.set<Real>(name) = in(name, 0.);
    name = "decay/a";          es.parameters.set<Real>(name) = in(name, 0.);
  }

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

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PIHNA");
  libmesh_assert_equal_to(system.n_vars(), 5);

  const System & struct_system =
    es.get_system<System>("uStructure");
  libmesh_assert_equal_to(struct_system.n_vars(), 2);

  FEType fe_type = system.variable_type(0);

  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
  //std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

  QGauss qrule (dim,   fe_type.default_quadrature_order());
  //QGauss qface (dim-1, fe_type.default_quadrature_order());

  fe->attach_quadrature_rule(&qrule);
  //fe_face->attach_quadrature_rule(&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  ///const std::vector<Real> & JxW_face = fe_face->get_JxW();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  //const std::vector<std::vector<Real>> & psi = fe_face->get_phi();

  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  //const std::vector<std::vector<RealGradient>> & dpsi = fe_face->get_dphi();

  //const std::vector<Point> & qface_points  = fe_face->get_xyz();
  //const std::vector<Point> & qface_normals = fe_face->get_normals();

  const Real DT_2 = es.parameters.get<Real>("time_step") / 2.0;

  const Real Lambda_k = es.parameters.get<Real>("cells_min_capacity");
  const Real Kappa_k = es.parameters.get<Real>("cells_max_capacity"),
             Kappa_a = es.parameters.get<Real>("cytokines_max_capacity");
  const Real ek = es.parameters.get<Real>("cells_max_capacity/exponent");
  const Real necrosis_c = es.parameters.get<Real>("necrosis/c") / Kappa_k,
             necrosis_h = es.parameters.get<Real>("necrosis/h") / Kappa_k,
             necrosis_v = es.parameters.get<Real>("necrosis/v") / Kappa_k;
  const Real diffuse_c_  = es.parameters.get<Real>("diffuse/c"),
             taxis_c_    = es.parameters.get<Real>("taxis/c"),
             diffuse_h_  = es.parameters.get<Real>("diffuse/h"),
             taxis_h_    = es.parameters.get<Real>("taxis/h"),
             produce_c   = es.parameters.get<Real>("produce/c"),
             switch_c2h  = es.parameters.get<Real>("switch/c/to/h"),
             switch_h2c  = es.parameters.get<Real>("switch/h/to/c"),
             switch_h2n  = es.parameters.get<Real>("switch/h/to/n");
  const Real diffuse_v_  = es.parameters.get<Real>("diffuse/v"),
             taxis_v_    = es.parameters.get<Real>("taxis/v"),
             produce_v   = es.parameters.get<Real>("produce/v");
  const Real secrete_a_c = es.parameters.get<Real>("secrete/a/from/c"),
             secrete_a_h = es.parameters.get<Real>("secrete/a/from/h"),
             uptake_a_v  = es.parameters.get<Real>("uptake/a/from/v"),
             decay_a     = es.parameters.get<Real>("decay/a");

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      system.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(5);
      for (unsigned int v=0; v<5; v++)
        system.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

      std::vector<std::vector<dof_id_type>> dof_indices_T_var(2);
      for (unsigned int l=0; l<2; l++)
        struct_system.get_dof_map().dof_indices(elem, dof_indices_T_var[l], l);

      const unsigned int n_dofs     = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      DenseMatrix<Number> Ke(n_dofs, n_dofs);
      DenseSubMatrix<Number> Ke_var[5][5] =
      {
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) }
      };
      for (unsigned int i=0; i<5; i++)
        for (unsigned int j=0; j<5; j++)
          Ke_var[i][j].reposition(i*n_var_dofs, j*n_var_dofs, n_var_dofs, n_var_dofs);

      DenseVector<Number> Fe(n_dofs);
      DenseSubVector<Number> Fe_var[5] =
      {
        DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)
      };
      for (unsigned int i=0; i<5; i++)
        Fe_var[i].reposition(i*n_var_dofs, n_var_dofs);

      fe->reinit(elem);

      struct_system.get_dof_map().dof_indices(elem, dof_indices_T_var[0], 0);
      //const Real HU = struct_system.solution->el(dof_indices_T_var[0][0]);
      struct_system.get_dof_map().dof_indices(elem, dof_indices_T_var[1], 1);
      //const Real RT = struct_system.solution->el(dof_indices_T_var[1][0]);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Number n_old(0.0), c_old(0.0), h_old(0.0), v_old(0.0), a_old(0.0);
          Gradient GRAD_c_old({0.0, 0.0, 0.0}), GRAD_h_old({0.0, 0.0, 0.0}), GRAD_v_old({0.0, 0.0, 0.0}), GRAD_a_old({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              n_old += phi[l][qp] * system.old_solution(dof_indices_var[0][l]);
              c_old += phi[l][qp] * system.old_solution(dof_indices_var[1][l]);
              h_old += phi[l][qp] * system.old_solution(dof_indices_var[2][l]);
              v_old += phi[l][qp] * system.old_solution(dof_indices_var[3][l]);
              a_old += phi[l][qp] * system.old_solution(dof_indices_var[4][l]);
              GRAD_c_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[1][l]));
              GRAD_h_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[2][l]));
              GRAD_v_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[3][l]));
              GRAD_a_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[4][l]));
            }

          const Real Tau = pow(1.0-apply_bounds(0.0, (n_old+c_old+h_old+v_old)/Kappa_k, 1.0), ek),
                     Tau__dn = -(ek*pow(1.0-apply_bounds(0.0, (n_old+c_old+h_old+v_old)/Kappa_k, 1.0), ek-1.0))/Kappa_k, Tau__dc = Tau__dn, Tau__dh = Tau__dn, Tau__dv = Tau__dn;

          const Real Ve = apply_bounds(0.0, v_old/(c_old+h_old+v_old), 1.0),
                     Ve__dc = -Ve/(c_old+h_old+v_old), Ve__dh = Ve__dc, Ve__dv = Ve/v_old+Ve__dc;

          const Real Ua = a_old/(a_old+Kappa_a),
                     Ua__da = 1.0/(a_old+Kappa_a)-Ua/(a_old+Kappa_a);

          const Real diffuse_c = (c_old>Lambda_k ? diffuse_c_ : 0.0),
                     taxis_c   = (c_old>Lambda_k ? taxis_c_   : 0.0),
                     diffuse_h = (h_old>Lambda_k ? diffuse_h_ : 0.0),
                     taxis_h   = (h_old>Lambda_k ? taxis_h_   : 0.0),
                     diffuse_v = (v_old>Lambda_k ? diffuse_v_ : 0.0),
                     taxis_v   = (v_old>Lambda_k ? taxis_v_   : 0.0);

          for (std::size_t i=0; i<n_var_dofs; i++)
            {
              // RHS contribution
              Fe_var[0](i) += JxW[qp]*(
                                        n_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               necrosis_c * c_old * n_old * phi[i][qp]
                                             + necrosis_h * h_old * n_old * phi[i][qp]
                                             + necrosis_v * v_old * n_old * phi[i][qp]
                                             + switch_h2n * (1.0-Ve) * h_old * phi[i][qp]
                                             )
                                      );
              // RHS contribution
              Fe_var[1](i) += JxW[qp]*(
                                        c_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               produce_c * Tau * c_old * phi[i][qp]
                                             - switch_c2h * (1.0-Ve) * c_old * phi[i][qp]
                                             + switch_h2c * Ve       * h_old * phi[i][qp]
                                             - necrosis_c * c_old * n_old * phi[i][qp]
                                             - diffuse_c * Tau * (GRAD_c_old * dphi[i][qp])
                                             - taxis_c * Tau * c_old * (GRAD_v_old * dphi[i][qp])
                                             )
                                      );
              // RHS contribution
              Fe_var[2](i) += JxW[qp]*(
                                        h_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               switch_c2h * (1.0-Ve) * c_old * phi[i][qp]
                                             - switch_h2c * Ve       * h_old * phi[i][qp]
                                             - necrosis_h * h_old * n_old * phi[i][qp]
                                             - diffuse_h * Tau * (GRAD_h_old * dphi[i][qp])
                                             - taxis_h * Tau * h_old * (GRAD_v_old * dphi[i][qp])
                                             - switch_h2n * (1.0-Ve) * h_old * phi[i][qp]
                                             )
                                      );
              // RHS contribution
              Fe_var[3](i) += JxW[qp]*(
                                        v_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               produce_v * Tau * Ua * v_old * phi[i][qp]
                                             - necrosis_v * v_old * n_old * phi[i][qp]
                                             - diffuse_v * Tau * (GRAD_v_old * dphi[i][qp])
                                             - taxis_v * Tau * v_old * (GRAD_a_old * dphi[i][qp])
                                             )
                                      );
              // RHS contribution
              Fe_var[4](i) += JxW[qp]*(
                                        a_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               secrete_a_c * c_old * phi[i][qp]
                                             + secrete_a_h * h_old * phi[i][qp]
                                             - uptake_a_v * v_old * a_old * phi[i][qp]
                                             - decay_a * a_old * phi[i][qp]
                                             )
                                      );

              for (std::size_t j=0; j<n_var_dofs; j++)
                {
                  // Matrix contribution
                  Ke_var[0][0](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        necrosis_c * c_old * phi[j][qp] * phi[i][qp]
                                                      + necrosis_h * h_old * phi[j][qp] * phi[i][qp]
                                                      + necrosis_v * v_old * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[0][1](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        necrosis_c * phi[j][qp] * n_old * phi[i][qp]
                                                      + switch_h2n * (-Ve__dc) * phi[j][qp] * h_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[0][2](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        necrosis_h * phi[j][qp] * n_old * phi[i][qp]
                                                      + switch_h2n * (-Ve__dh) * phi[j][qp] * h_old * phi[i][qp]
                                                      + switch_h2n * (1.0-Ve) * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[0][3](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        necrosis_v * phi[j][qp] * n_old * phi[i][qp]
                                                      + switch_h2n * (-Ve__dv) * phi[j][qp] * h_old * phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[1][0](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_c * Tau__dn * phi[j][qp] * c_old * phi[i][qp]
                                                      - necrosis_c * c_old * phi[j][qp] * phi[i][qp]
                                                      - diffuse_c * Tau__dn * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - taxis_c * Tau__dn * phi[j][qp] * c_old * (GRAD_v_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[1][1](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        produce_c * Tau * phi[j][qp] * phi[i][qp]
                                                      + produce_c * Tau__dc * phi[j][qp] * c_old * phi[i][qp]
                                                      - switch_c2h * (1.0-Ve) * phi[j][qp] * phi[i][qp]
                                                      - switch_c2h * (-Ve__dc) * phi[j][qp] * c_old * phi[i][qp]
                                                      + switch_h2c * Ve__dc    * phi[j][qp] * h_old * phi[i][qp]
                                                      - necrosis_c * phi[j][qp] * n_old * phi[i][qp]
                                                      - diffuse_c * Tau__dc * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - diffuse_c * Tau * (dphi[j][qp] * dphi[i][qp])
                                                      - taxis_c * Tau__dc * phi[j][qp] * c_old * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_c * Tau * phi[j][qp] * (GRAD_v_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[1][2](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_c * Tau__dh * phi[j][qp] * c_old * phi[i][qp]
                                                      - switch_c2h * (-Ve__dh) * phi[j][qp] * c_old * phi[i][qp]
                                                      + switch_h2c * Ve__dh    * phi[j][qp] * h_old * phi[i][qp]
                                                      + switch_h2c * Ve        * phi[j][qp] * phi[i][qp]
                                                      - diffuse_c * Tau__dh * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - taxis_c * Tau__dh * phi[j][qp] * c_old * (GRAD_v_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[1][3](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_c * Tau__dv * phi[j][qp] * c_old * phi[i][qp]
                                                      - switch_c2h * (-Ve__dv) * phi[j][qp] * c_old * phi[i][qp]
                                                      + switch_h2c * Ve__dv    * phi[j][qp] * h_old * phi[i][qp]
                                                      - diffuse_c * Tau__dv * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - taxis_c * Tau__dv * phi[j][qp] * c_old * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_c * Tau * c_old * (dphi[j][qp] * dphi[i][qp])
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[2][0](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                      - necrosis_h * h_old * phi[j][qp] * phi[i][qp]
                                                      - diffuse_h * Tau__dn * phi[j][qp] * (GRAD_h_old * dphi[i][qp])
                                                      - taxis_h * Tau__dn * phi[j][qp] * h_old * (GRAD_v_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[2][1](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        switch_c2h * (1.0-Ve) * phi[j][qp] * phi[i][qp]
                                                      + switch_c2h * (-Ve__dc) * phi[j][qp] * c_old * phi[i][qp]
                                                      - switch_h2c * Ve__dc    * phi[j][qp] * h_old * phi[i][qp]
                                                      - diffuse_h * Tau__dc * phi[j][qp] * (GRAD_h_old * dphi[i][qp])
                                                      - taxis_h * Tau__dc * phi[j][qp] * h_old * (GRAD_v_old * dphi[i][qp])
                                                      - switch_h2n * (-Ve__dc) * phi[j][qp] * h_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][2](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        switch_c2h * (-Ve__dh) * phi[j][qp] * c_old * phi[i][qp]
                                                      - switch_h2c * Ve__dh    * phi[j][qp] * h_old * phi[i][qp]
                                                      - switch_h2c * Ve       * phi[j][qp] * phi[i][qp]
                                                      - necrosis_h * phi[j][qp] * n_old * phi[i][qp]
                                                      - diffuse_h * Tau__dh * phi[j][qp] * (GRAD_h_old * dphi[i][qp])
                                                      - diffuse_h * Tau * (dphi[j][qp] * dphi[i][qp])
                                                      - taxis_h * Tau__dh * phi[j][qp] * h_old * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_h * Tau * phi[j][qp] * (GRAD_v_old * dphi[i][qp])
                                                      - switch_h2n * (-Ve__dh) * phi[j][qp] * h_old * phi[i][qp]
                                                      - switch_h2n * (1.0-Ve) * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][3](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        switch_c2h * (-Ve__dv) * phi[j][qp] * c_old * phi[i][qp]
                                                      - switch_h2c * Ve__dv    * phi[j][qp] * h_old * phi[i][qp]
                                                      - diffuse_h * Tau__dv * phi[j][qp] * (GRAD_h_old * dphi[i][qp])
                                                      - taxis_h * Tau__dv * phi[j][qp] * h_old * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_h * Tau * h_old * (dphi[j][qp] * dphi[i][qp])
                                                      - switch_h2n * (-Ve__dv) * phi[j][qp] * h_old * phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[3][0](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_v * Tau__dn * phi[j][qp] * Ua * v_old * phi[i][qp]
                                                      - necrosis_v * v_old * phi[j][qp] * phi[i][qp]
                                                      - diffuse_v * Tau__dn * phi[j][qp] * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_v * Tau__dn * phi[j][qp] * v_old * (GRAD_a_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[3][1](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_v * Tau__dc * phi[j][qp] * Ua * v_old * phi[i][qp]
                                                      - diffuse_v * Tau__dc * phi[j][qp] * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_v * Tau__dc * phi[j][qp] * v_old * (GRAD_a_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[3][2](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_v * Tau__dh * phi[j][qp] * Ua * v_old * phi[i][qp]
                                                      - diffuse_v * Tau__dh * phi[j][qp] * (GRAD_v_old * dphi[i][qp])
                                                      - taxis_v * Tau__dh * phi[j][qp] * v_old * (GRAD_a_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[3][3](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        produce_v * Tau__dv * phi[j][qp] * Ua * v_old * phi[i][qp]
                                                      - necrosis_v * phi[j][qp] * n_old * phi[i][qp]
                                                      - diffuse_v * Tau__dv * phi[j][qp] * (GRAD_v_old * dphi[i][qp])
                                                      - diffuse_v * Tau * (dphi[j][qp] * dphi[i][qp])
                                                      - taxis_v * Tau__dv * phi[j][qp] * v_old * (GRAD_a_old * dphi[i][qp])
                                                      - taxis_v * Tau * phi[j][qp] * (GRAD_a_old * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[3][4](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        produce_v * Tau * Ua__da * phi[j][qp] * v_old * phi[i][qp]
                                                      - taxis_v * Tau * v_old * (dphi[j][qp] * dphi[i][qp])
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[4][1](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        secrete_a_c * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[4][2](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        secrete_a_h * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[4][3](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                      - uptake_a_v * phi[j][qp] * a_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[4][4](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                      - uptake_a_v * v_old * phi[j][qp] * phi[i][qp]
                                                      - decay_a * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                }
            }
        }

      system.get_dof_map().constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      //
      system.get_system_matrix().add_matrix(Ke, dof_indices);
      system.rhs->add_vector(Fe, dof_indices);
    }
  // ...done
}

void check_solution (EquationSystems & es)
{
  const MeshBase& mesh = es.get_mesh();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PIHNA");
  libmesh_assert_equal_to(system.n_vars(), 5);

  std::vector<Number> soln;
  system.update_global_solution(soln);

  for (const auto & node : mesh.node_ptr_range())
    {
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

      Real n_, c_, h_, v_, a_;
      n_ = soln[idof[0]]; if (n_<0.0) n_ = 0.0;
      c_ = soln[idof[1]]; if (c_<0.0) c_ = 0.0;
      h_ = soln[idof[2]]; if (h_<0.0) h_ = 0.0;
      v_ = soln[idof[3]]; if (v_<0.0) v_ = 0.0;
      a_ = soln[idof[4]]; if (a_<0.0) a_ = 0.0;

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
