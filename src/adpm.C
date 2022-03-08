#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/equation_systems.h"
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
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

void input (GetPot & , EquationSystems & );
void initial_adpm (EquationSystems & , const std::string & );
void assemble_adpm (EquationSystems & , const std::string & );
void initial_tracts (EquationSystems & , const std::string & );

extern PerfLog plog;

void adpm (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);
  GetPot in("input.dat");

  input(in, es);

  TransientLinearImplicitSystem & adpm =
    es.add_system<TransientLinearImplicitSystem>("ADPM");
  adpm.add_variable("PrP", FIRST, LAGRANGE);
  adpm.add_variable("A_b", FIRST, LAGRANGE);
  adpm.add_variable("Tau", FIRST, LAGRANGE);
  adpm.attach_assemble_function(assemble_adpm);
  adpm.attach_init_function(initial_adpm);

  ExplicitSystem & tracts =
    es.add_system<ExplicitSystem>("Tracts");
  tracts.add_variable("TractX", CONSTANT, MONOMIAL);
  tracts.add_variable("TractY", CONSTANT, MONOMIAL);
  tracts.add_variable("TractZ", CONSTANT, MONOMIAL);
  tracts.attach_init_function(initial_tracts);

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
      adpm.time = es.parameters.get<Real>("time");
      const Real current_time = adpm.time;

      libMesh::out << " Solving time increment: " << t
                   << " (time=" << current_time <<  ") ..." << std::endl;

      // copy the previously-current solution into the old solution
      *(adpm.old_local_solution) = *(adpm.current_local_solution);
      // now solve the AD progression model
      adpm.solve();

      if (0 == t%output_step)
        ex2.write_timestep(ex2_filename, es, t, current_time);
    }

  // ...done
}

void input (GetPot & in, EquationSystems & es)
{
  std::string name;

  name = "input_GMSH";
  es.parameters.set<std::string>(name) = in(name, "input.msh");
  name = "output_GMSH";
  es.parameters.set<std::string>(name) = in(name, "output.msh");
  name = "input_nodal";
  es.parameters.set<std::string>(name) = in(name, "input.nodal");
  name = "input_elemental";
  es.parameters.set<std::string>(name) = in(name, "input.elemental");
  name = "output_EXODUS";
  es.parameters.set<std::string>(name) = in(name, "output.ex2");

  es.parameters.set<RealVectorValue>("velocity") = RealVectorValue(0., 0., 0.);

  es.parameters.set<Real>("time") = 0.0;

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = in(name, 1);
  name = "output_step";
  es.parameters.set<int>(name) = in(name, 1);

  // parameters for the species: PrP
  {
    name = "decay/PrP";             es.parameters.set<Real>(name) = in(name, 0.);
    name = "decay/PrP/pulse/0";     es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "decay/PrP/pulse/1";     es.parameters.set<Real>(name) = in(name,+1.0e+20);
  }

  // parameters for the species: A_b
  {
    name = "diffuse/A_b";           es.parameters.set<Real>(name) = in(name, 0.);
    name = "diffuse/A_b/pulse/0";   es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "diffuse/A_b/pulse/1";   es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "taxis/A_b";             es.parameters.set<Real>(name) = in(name, 0.);
    name = "taxis/A_b/pulse/0";     es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "taxis/A_b/pulse/1";     es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "produce/A_b";           es.parameters.set<Real>(name) = in(name, 0.);
    name = "produce/A_b/sigmoid/0"; es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "produce/A_b/sigmoid/1"; es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "transform/A_b";         es.parameters.set<Real>(name) = in(name, 0.);
    name = "transform/A_b/pulse/0"; es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "transform/A_b/pulse/1"; es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "decay/A_b";             es.parameters.set<Real>(name) = in(name, 0.);
    name = "decay/A_b/pulse/0";     es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "decay/A_b/pulse/1";     es.parameters.set<Real>(name) = in(name,+1.0e+20);
  }

  // parameters for the species: Tau
  {
    name = "diffuse/Tau";           es.parameters.set<Real>(name) = in(name, 0.);
    name = "diffuse/Tau/pulse/0";   es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "diffuse/Tau/pulse/1";   es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "taxis/Tau";             es.parameters.set<Real>(name) = in(name, 0.);
    name = "taxis/Tau/pulse/0";     es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "taxis/Tau/pulse/1";     es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "produce/Tau";           es.parameters.set<Real>(name) = in(name, 0.);
    name = "produce/Tau/sigmoid/0"; es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "produce/Tau/sigmoid/1"; es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "transform/Tau";         es.parameters.set<Real>(name) = in(name, 0.);
    name = "transform/Tau/pulse/0"; es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "transform/Tau/pulse/1"; es.parameters.set<Real>(name) = in(name,+1.0e+20);
    name = "decay/Tau";             es.parameters.set<Real>(name) = in(name, 0.);
    name = "decay/Tau/pulse/0";     es.parameters.set<Real>(name) = in(name,-1.0e-20);
    name = "decay/Tau/pulse/1";     es.parameters.set<Real>(name) = in(name,+1.0e+20);
  }

  // ...done
}

void initial_tracts (EquationSystems & es,
                     const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "Tracts");

  const MeshBase& mesh = es.get_mesh();

  ExplicitSystem & system =
    es.get_system<ExplicitSystem>("Tracts");
  libmesh_assert_equal_to(system.n_vars(), 3);

  std::ifstream fin(es.parameters.get<std::string>("input_elemental"));

  for (const auto & elem : mesh.active_element_ptr_range())
    {
      Real vec[3];
      fin >> vec[0] >> vec[1] >> vec[2];

      std::vector<std::vector<dof_id_type>> dof_indices_T_var(3);

      for (unsigned int l=0; l<3; l++)
        {
          system.get_dof_map().dof_indices(elem, dof_indices_T_var[l], l);
          system.solution->set(dof_indices_T_var[l][0], vec[l]);
        }
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

inline // Pi or rectangular function
Real Pi_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  0.0;
  else if (C < c1) return  cM;
  else             return  0.0;
}

inline // step-decay function
Real SD_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  cM;
  else if (C < c1) return  cM*(c1-C)/(c1-c0);
  else             return  0.0;
}

inline // step-growth function
Real SG_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & c0 = p_[1];
  const Real & c1 = p_[2];
  if      (C < c0) return  cM;
  else if (C < c1) return  cM*(C-c0)/(c1-c0);
  else             return  0.0;
}

inline // Boltzmann (sigmoidal) increase function
Real Bsi_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return G/(1.0+G);
}
inline // Boltzmann (sigmoidal) increase function - derivative
Real deriv_Bsi_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return G/(dC*(1.0+G)*(1.0+G));
}

inline // Boltzmann (sigmoidal) dencrease function
Real Bsd_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return 1.0/(1.0+G);
}
inline // Boltzmann (sigmoidal) dencrease function - derivative
Real deriv_Bsd_(const Real & C, const Real * p_)
{
  const Real & cM = p_[0];
  if (0.0==cM) return 0.0;
  const Real & C0 = p_[1];
  const Real & dC = p_[2];
  const Real G = exp((C-C0)/dC);
  return -G/(dC*(1.0+G)*(1.0+G));
}

void initial_adpm (EquationSystems & es,
                   const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "ADPM");

  const MeshBase& mesh = es.get_mesh();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("ADPM");
  libmesh_assert_equal_to(system.n_vars(), 3);

  es.parameters.set<Real> ("time") =
  system.time = 0.0;

  std::ifstream fin(es.parameters.get<std::string>("input_nodal"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real PrP, A_b, Tau;
      fin >> PrP >> A_b >> Tau;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );

      system.solution->set(idof[0], PrP);
      system.solution->set(idof[1], A_b);
      system.solution->set(idof[2], Tau);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void assemble_adpm (EquationSystems & es,
                    const std::string & system_name)
{
  libmesh_ignore(es, system_name);
  libmesh_assert_equal_to(system_name, "ADPM");
  libmesh_assert_equal_to(system.n_vars(), 3);

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("ADPM");

  FEType fe_type = system.variable_type(0);

  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
  std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

  QGauss qrule (dim,   fe_type.default_quadrature_order());
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  fe->attach_quadrature_rule(&qrule);
  fe_face->attach_quadrature_rule(&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<Real> & JxW_face = fe_face->get_JxW();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<Real>> & psi = fe_face->get_phi();

  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  const std::vector<std::vector<RealGradient>> & dpsi = fe_face->get_dphi();

  //const std::vector<Point> & qface_points  = fe_face->get_xyz();
  const std::vector<Point> & qface_normals = fe_face->get_normals();

  const System & tracts_system = es.get_system<System>("Tracts");

  const Real DT_2 = es.parameters.get<Real>("time_step") / 2.0;

  //const RealVectorValue velocity = es.parameters.get<RealVectorValue>("velocity");

  const Real decay_PrP[] = { es.parameters.get<Real>("decay/PrP")         ,
                             es.parameters.get<Real>("decay/PrP/pulse/0") ,
                             es.parameters.get<Real>("decay/PrP/pulse/1") };
  const Real diffuse_A_b[] = { es.parameters.get<Real>("diffuse/A_b")         ,
                               es.parameters.get<Real>("diffuse/A_b/pulse/0") ,
                               es.parameters.get<Real>("diffuse/A_b/pulse/1") };
  const Real taxis_A_b[]   = { es.parameters.get<Real>("taxis/A_b")         ,
                               es.parameters.get<Real>("taxis/A_b/pulse/0") ,
                               es.parameters.get<Real>("taxis/A_b/pulse/1") };
  const Real produce_A_b[] = { es.parameters.get<Real>("produce/A_b")           ,
                               es.parameters.get<Real>("produce/A_b/sigmoid/0") ,
                               es.parameters.get<Real>("produce/A_b/sigmoid/1") };
  const Real transform_A_b[] = { es.parameters.get<Real>("transform/A_b")         ,
                                 es.parameters.get<Real>("transform/A_b/pulse/0") ,
                                 es.parameters.get<Real>("transform/A_b/pulse/1") };
  const Real decay_A_b[] = { es.parameters.get<Real>("decay/A_b")         ,
                             es.parameters.get<Real>("decay/A_b/pulse/0") ,
                             es.parameters.get<Real>("decay/A_b/pulse/1") };
  const Real diffuse_Tau[] = { es.parameters.get<Real>("diffuse/Tau")         ,
                               es.parameters.get<Real>("diffuse/Tau/pulse/0") ,
                               es.parameters.get<Real>("diffuse/Tau/pulse/1") };
  const Real taxis_Tau[]   = { es.parameters.get<Real>("taxis/Tau")         ,
                               es.parameters.get<Real>("taxis/Tau/pulse/0") ,
                               es.parameters.get<Real>("taxis/Tau/pulse/1") };
  const Real produce_Tau[] = { es.parameters.get<Real>("produce/Tau")           ,
                               es.parameters.get<Real>("produce/Tau/sigmoid/0") ,
                               es.parameters.get<Real>("produce/Tau/sigmoid/1") };
  const Real transform_Tau[] = { es.parameters.get<Real>("transform/Tau")         ,
                                 es.parameters.get<Real>("transform/Tau/pulse/0") ,
                                 es.parameters.get<Real>("transform/Tau/pulse/1") };
  const Real decay_Tau[] = { es.parameters.get<Real>("decay/Tau")         ,
                             es.parameters.get<Real>("decay/Tau/pulse/0") ,
                             es.parameters.get<Real>("decay/Tau/pulse/1") };

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      system.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(3);
      for (unsigned int v=0; v<3; v++)
        system.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

      std::vector<std::vector<dof_id_type>> dof_indices_T_var(3);
      for (unsigned int l=0; l<3; l++)
        system.get_dof_map().dof_indices(elem, dof_indices_T_var[l], l);

      const unsigned int n_dofs     = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      DenseMatrix<Number> Ke(n_dofs, n_dofs);
      DenseSubMatrix<Number> Ke_var[3][3] =
      {
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
        { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) }
      };
      for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
          Ke_var[i][j].reposition(i*n_var_dofs, j*n_var_dofs, n_var_dofs, n_var_dofs);

      DenseVector<Number> Fe(n_dofs);
      DenseSubVector<Number> Fe_var[3] =
      {
        DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)
      };
      for (unsigned int i=0; i<3; i++)
        Fe_var[i].reposition(i*n_var_dofs, n_var_dofs);

      fe->reinit(elem);

      Point tracts;
      for (unsigned int l=0; l<3; l++)
        {
          tracts_system.get_dof_map().dof_indices(elem, dof_indices_T_var[l], l);
          tracts(l) = tracts_system.solution->el(dof_indices_T_var[l][0]);
        }

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          VectorValue<Number>        u_old(0.0, 0.0, 0.0);
          VectorValue<Gradient> GRAD_u_old({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              u_old(0) += phi[l][qp] * system.old_solution(dof_indices_var[0][l]);
              u_old(1) += phi[l][qp] * system.old_solution(dof_indices_var[1][l]);
              u_old(2) += phi[l][qp] * system.old_solution(dof_indices_var[2][l]);
              GRAD_u_old(1).add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[1][l]));
              GRAD_u_old(2).add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[2][l]));
            }

          for (std::size_t i=0; i<n_var_dofs; i++)
            {
              // RHS contribution
              Fe_var[0](i) += JxW[qp]*(
                                        // capacity (mass) term
                                        phi[i][qp]*u_old(0)
                                        // convection, diffusion, source/sink terms
                                      + DT_2*(
                                             - Pi_(u_old(1),transform_A_b)*u_old(0)*phi[i][qp]
                                             - Pi_(u_old(2),transform_Tau)*u_old(0)*phi[i][qp]
                                             - Pi_(u_old(0),decay_PrP)*u_old(0)*phi[i][qp]
                                             //- (GRAD_u_old(0)*velocity)*phi[i][qp]
                                             )
                                      );
              Fe_var[1](i) += JxW[qp]*(
                                        // capacity (mass) term
                                        phi[i][qp]*u_old(1)
                                        // convection, diffusion, source/sink terms
                                      + DT_2*(
                                               SD_(u_old(1),produce_A_b)*u_old(1)*phi[i][qp]
                                             + Pi_(u_old(1),transform_A_b)*u_old(0)*phi[i][qp]
                                             - Pi_(u_old(1),decay_A_b)*u_old(1)*phi[i][qp]
                                             - Pi_(u_old(1),diffuse_A_b)*(GRAD_u_old(1)*dphi[i][qp])
                                             - Pi_(u_old(1),taxis_A_b)*((GRAD_u_old(1)*tracts)*tracts*dphi[i][qp])
                                             //- (GRAD_u_old(1)*velocity)*phi[i][qp]
                                             )
                                      );
              Fe_var[2](i) += JxW[qp]*(
                                        // capacity (mass) term
                                        phi[i][qp]*u_old(2)
                                        // convection, diffusion, source/sink terms
                                      + DT_2*(
                                               SD_(u_old(2),produce_Tau)*u_old(2)*phi[i][qp]
                                             + Pi_(u_old(2),transform_Tau)*u_old(0)*phi[i][qp]
                                             - Pi_(u_old(2),decay_Tau)*u_old(2)*phi[i][qp]
                                             - Pi_(u_old(2),diffuse_Tau)*(GRAD_u_old(2)*dphi[i][qp])
                                             - Pi_(u_old(2),taxis_Tau)*((GRAD_u_old(2)*tracts)*tracts*dphi[i][qp])
                                             //- (GRAD_u_old(2)*velocity)*phi[i][qp]
                                             )
                                      );

              for (std::size_t j=0; j<n_var_dofs; j++)
                {
                  // Matrix contribution
                  Ke_var[0][0](i,j) += JxW[qp]*(
                                                 // capacity (mass) term
                                                 phi[i][qp]*phi[j][qp]
                                                 // convection, diffusion, source/sink terms
                                               - DT_2*(
                                                      - Pi_(u_old(1),transform_A_b)*phi[i][qp]*phi[j][qp]
                                                      - Pi_(u_old(2),transform_Tau)*phi[i][qp]*phi[j][qp]
                                                      - Pi_(u_old(0),decay_PrP)*phi[i][qp]*phi[j][qp]
                                                      //- (dphi[j][qp]*velocity)*phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[1][0](i,j) += JxW[qp]*(
                                                 // convection, diffusion, source/sink terms
                                               - DT_2*(
                                                        Pi_(u_old(1),transform_A_b)*phi[i][qp]*phi[j][qp]
                                                      )
                                               );
                  Ke_var[1][1](i,j) += JxW[qp]*(
                                                 // capacity (mass) term
                                                 phi[i][qp]*phi[j][qp]
                                                 // convection, diffusion, source/sink terms
                                               - DT_2*(
                                                        SD_(u_old(1),produce_A_b)*phi[i][qp]*phi[j][qp]
                                                      - Pi_(u_old(1),decay_A_b)*phi[i][qp]*phi[j][qp]
                                                      - Pi_(u_old(1),diffuse_A_b)*(dphi[i][qp]*dphi[j][qp])
                                                      - Pi_(u_old(1),taxis_A_b)*((dphi[j][qp]*tracts)*tracts*dphi[i][qp])
                                                      //- (dphi[j][qp]*velocity)*phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[2][0](i,j) += JxW[qp]*(
                                                 // convection, diffusion, source/sink terms
                                               - DT_2*(
                                                        Pi_(u_old(2),transform_Tau)*phi[i][qp]*phi[j][qp]
                                                      )
                                               );
                  Ke_var[2][2](i,j) += JxW[qp]*(
                                                 // capacity (mass) term
                                                 phi[i][qp]*phi[j][qp]
                                                 // convection, diffusion, source/sink terms
                                               - DT_2*(
                                                        SD_(u_old(2),produce_Tau)*phi[i][qp]*phi[j][qp]
                                                      - Pi_(u_old(2),decay_Tau)*phi[i][qp]*phi[j][qp]
                                                      - Pi_(u_old(2),diffuse_Tau)*(dphi[i][qp]*dphi[j][qp])
                                                      - Pi_(u_old(2),taxis_Tau)*((dphi[j][qp]*tracts)*tracts*dphi[i][qp])
                                                      //- (dphi[j][qp]*velocity)*phi[i][qp]
                                                      )
                                               );
                }
            }
        }

      if ( 0 )
      {
        // penalty value...
        const Real penalty = 1.0e+10;

        for (auto s : elem->side_index_range())
          {
            if (elem->neighbor_ptr(s) != nullptr) continue;

            fe_face->reinit(elem, s);

            for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {

                VectorValue<Number>        u_old(0.0, 0.0, 0.0);
                VectorValue<Gradient> GRAD_u_old({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
                for (std::size_t l=0; l<n_var_dofs; l++)
                  {
                    u_old(1) += psi[l][qp] * system.old_solution(dof_indices_var[1][l]);
                    u_old(2) += psi[l][qp] * system.old_solution(dof_indices_var[2][l]);
                    GRAD_u_old(1).add_scaled(dpsi[l][qp], system.old_solution(dof_indices_var[1][l]));
                    GRAD_u_old(2).add_scaled(dpsi[l][qp], system.old_solution(dof_indices_var[2][l]));
                  }

                VectorValue<Number> field(0.0, 0.0, 0.0);

                field(1) = qface_normals[qp] * ( //velocity*u_old(1) +
                                                 Pi_(u_old(1),diffuse_A_b)*GRAD_u_old(1));
                field(2) = qface_normals[qp] * ( //velocity*u_old(2) +
                                                 Pi_(u_old(2),diffuse_Tau)*GRAD_u_old(2));

                for (std::size_t i=0; i<psi.size(); i++)
                  {
                    // RHS contribution
                    Fe_var[1](i) += JxW_face[qp] * (psi[i][qp]*field(1))
                                  * penalty;
                    Fe_var[2](i) += JxW_face[qp] * (psi[i][qp]*field(2))
                                  * penalty;
                    // Matrix contribution
                    for (std::size_t j=0; j<psi.size(); j++)
                      {
                        Ke_var[1][1](i,j) += JxW_face[qp] * (psi[i][qp]*psi[j][qp])
                                           * penalty;
                        Ke_var[2][2](i,j) += JxW_face[qp] * (psi[i][qp]*psi[j][qp])
                                           * penalty;
                      }
                  }
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
