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

static void input (const std::string & , EquationSystems & );
static void initial_radiotherapy (EquationSystems & , const std::string & );
static void initial_ripf (EquationSystems & , const std::string & );
static void assemble_ripf (EquationSystems & , const std::string & );
static void check_solution (EquationSystems & );

extern PerfLog plog;

void ripf (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);

  input("input.dat", es);

  TransientLinearImplicitSystem & model =
    es.add_system<TransientLinearImplicitSystem>("RIPF");
  model.add_variable("HU", FIRST, LAGRANGE);
  model.add_variable("cc", FIRST, LAGRANGE);
  model.add_variable("fb", FIRST, LAGRANGE);
  model.attach_assemble_function(assemble_ripf);
  model.attach_init_function(initial_ripf);

  ExplicitSystem & radiotherapy =
    es.add_system<ExplicitSystem>("RT");
  radiotherapy.add_variable("RT_dose/broad", FIRST, LAGRANGE);
  radiotherapy.add_variable("RT_dose/focus", FIRST, LAGRANGE);
  radiotherapy.add_variable("RT_dose/total", FIRST, LAGRANGE);
  radiotherapy.attach_init_function(initial_radiotherapy);

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
      *(model.older_local_solution) = *(model.old_local_solution);
      *(model.old_local_solution) = *(model.current_local_solution);
      // now solve the AD progression model
      //++++++++++++++++++++++++++++++++++++++model.solve();

      check_solution(es);

      if (0 == t%output_step)
        ex2.write_timestep(ex2_filename, es, t, current_time);
    }

  // ...done
}

void input (const std::string & file_name, EquationSystems & es)
{
  GetPot in(file_name);

  std::string name;

  // create a time-stamped directory to store in simulation results
  const std::string DIR = date_time_to_string(date_now(), "%Y%m%d_%H%M%S") + "/";
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
  name = "input_nodal_RT";
  es.parameters.set<std::string>(name) = in(name, "input.nodal~RT");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());
  //
  name = "output_EXODUS";
  es.parameters.set<std::string>(name) = DIR + in(name, "output.ex2");

  es.parameters.set<RealVectorValue>("velocity") = RealVectorValue(0., 0., 0.);

  es.parameters.set<Real>("time") = 0.0;

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = in(name, 1);
  name = "output_step";
  es.parameters.set<int>(name) = in(name, 1);

  // general parameters including for radiotherapy (RT)
  {
    name = "RT_dose/broad/fractions"; es.parameters.set<int>(name) = in(name, 1);
    name = "RT_dose/focus/fractions"; es.parameters.set<int>(name) = in(name, 1);
    //
    name = "volume_fraction/stroma";     es.parameters.set<Real>(name) = in(name, 0.);
    name = "volume_fraction/parenchyma"; es.parameters.set<Real>(name) = in(name, 0.);
    name = "volume_fraction/exponent";   es.parameters.set<Real>(name) = in(name, 1.);
    name = "volume_fraction/min_vacant"; es.parameters.set<Real>(name) = in(name, 1.e-12);
    const Real VF_ = es.parameters.get<Real>(name);
    name = "volume_fraction/max_vacant"; es.parameters.set<Real>(name) = in(name, 1.-VF_);
    //
    name = "HU/min"; es.parameters.set<Real>(name) = in(name, -1000.);
    name = "HU/max"; es.parameters.set<Real>(name) = in(name, +1000.);
  }

  // parameters for the species: cc
  {
    name = "cc/kappa";      es.parameters.set<Real>(name) = in(name, 0.);
    name = "cc/delta";      es.parameters.set<Real>(name) = in(name, 0.);
    name = "cc/delta/RT/a"; es.parameters.set<Real>(name) = in(name, 0.);
    name = "cc/delta/RT/b"; es.parameters.set<Real>(name) = in(name, 0.);
  }

  // parameters for the species: fb
  {
    name = "fb/lambda";      es.parameters.set<Real>(name) = in(name, 0.);
    name = "fb/lambda/RT/r"; es.parameters.set<Real>(name) = in(name, 1.);
    name = "fb/diffusion";   es.parameters.set<Real>(name) = in(name, 0.);
    name = "fb/haptotaxis";  es.parameters.set<Real>(name) = in(name, 0.);
  }

  // ...done
}

void initial_radiotherapy (EquationSystems & es,
                           const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "RT");

  const MeshBase& mesh = es.get_mesh();

  ExplicitSystem & system =
    es.get_system<ExplicitSystem>("RT");
  libmesh_assert_equal_to(system.n_vars(), 3);

  std::ifstream fin(es.parameters.get<std::string>("input_nodal_RT"));

  const Real RT_broad_frac = es.parameters.get<int>("RT_dose/broad/fractions"),
             RT_focus_frac = es.parameters.get<int>("RT_dose/focus/fractions"),
             RT_total_frac = RT_broad_frac + RT_focus_frac;
  // simulation time is expressed in "days"
  const int day = std::floor( 0.0 );

  for (const auto & node : mesh.node_ptr_range())
    {
      Real RT_broad_, RT_focus_;
      fin >> RT_broad_ >> RT_focus_;
      Real RT_total_ = 0.0; // radiation therapy (RT) dose per fraction
      if      ( day < RT_broad_frac ) RT_total_ = RT_broad_ / RT_broad_frac * (day+1);
      else if ( day < RT_total_frac ) RT_total_ = RT_focus_ / RT_focus_frac * ((day+1)-RT_broad_frac) + RT_broad_;
      else                            RT_total_ = RT_broad_ + RT_focus_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );

      system.solution->set(idof[0], RT_broad_);
      system.solution->set(idof[1], RT_focus_);
      system.solution->set(idof[2], RT_total_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void initial_ripf (EquationSystems & es,
                   const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "RIPF");

  const MeshBase& mesh = es.get_mesh();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("RIPF");
  libmesh_assert_equal_to(system.n_vars(), 3);

  es.parameters.set<Real> ("time") =
  system.time = 0.0;

  std::ifstream fin(es.parameters.get<std::string>("input_nodal"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real HU_, cc_, fb_;
      fin >> HU_ >> cc_ >> fb_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );

      system.solution->set(idof[0], HU_);
      system.solution->set(idof[1], cc_);
      system.solution->set(idof[2], fb_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void assemble_ripf (EquationSystems & es,
                    const std::string & system_name)
{
  libmesh_ignore(es, system_name);
  libmesh_assert_equal_to(system_name, "RIPF");

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("RIPF");
  libmesh_assert_equal_to(system.n_vars(), 3);

  const ExplicitSystem & RT_system =
    es.get_system<ExplicitSystem>("RT");
  libmesh_assert_equal_to(RT_system.n_vars(), 3);

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

  const Real DT_2 = es.parameters.get<Real>("time_step") / 2.0;
  const Real DT_R = 1.0 / es.parameters.get<Real>("time_step");

  const Real VolFr_stroma     = es.parameters.get<Real>("volume_fraction/stroma"),
             VolFr_parenchyma = es.parameters.get<Real>("volume_fraction/parenchyma"),
             VolFr_exponent   = es.parameters.get<Real>("volume_fraction/exponent"),
             VolFr_min_vacant = es.parameters.get<Real>("volume_fraction/min_vacant"),
             VolFr_max_vacant = es.parameters.get<Real>("volume_fraction/max_vacant");

  const Real kappa      = es.parameters.get<Real>("cc/kappa"),
             delta      = es.parameters.get<Real>("cc/delta"),
             delta_RT_a = es.parameters.get<Real>("cc/delta/RT/a"),
             delta_RT_b = es.parameters.get<Real>("cc/delta/RT/b");
  const Real lambda      = es.parameters.get<Real>("fb/lambda"),
             lambda_RT_r = es.parameters.get<Real>("fb/lambda/RT/r"),
             diffusion   = es.parameters.get<Real>("fb/diffusion"),
             haptotaxis  = es.parameters.get<Real>("fb/haptotaxis");

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      system.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(3);
      for (unsigned int v=0; v<3; v++)
        system.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

      std::vector<std::vector<dof_id_type>> dof_indices_RT_var(3);
      for (unsigned int l=0; l<3; l++)
        RT_system.get_dof_map().dof_indices(elem, dof_indices_RT_var[l], l);

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

      const Real RT = RT_system.solution->el(dof_indices_RT_var[2][0]);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Number HU_old(0.0), cc_old(0.0), fb_old(0.0);
          Gradient GRAD_HU_old({0.0, 0.0, 0.0}), GRAD_fb_old({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              HU_old += phi[l][qp] * system.old_solution(dof_indices_var[0][l]);
              cc_old += phi[l][qp] * system.old_solution(dof_indices_var[1][l]);
              fb_old += phi[l][qp] * system.old_solution(dof_indices_var[2][l]);
              GRAD_HU_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[0][l]));
              GRAD_fb_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[2][l]));
            }
          Number cc_older(0.0), fb_older(0.0);
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              cc_older += phi[l][qp] * system.older_solution(dof_indices_var[1][l]);
              fb_older += phi[l][qp] * system.older_solution(dof_indices_var[2][l]);
            }
          Number cc__dtime((cc_old-cc_older)*DT_R), fb__dtime((fb_old-fb_older)*DT_R);

          const Real VolFr_cells = cc_old + fb_old;
          const Real VolFr_TOTAL = VolFr_stroma + VolFr_parenchyma + VolFr_cells;

          Real Tau = 0.0;
          Real Tau__dcc = 0.0,
               Tau__dfb = 0.0;
          if (1.0 > VolFr_TOTAL)
            {
              Tau = pow(1.0-VolFr_TOTAL, VolFr_exponent);
              Tau__dcc =
              Tau__dfb = -VolFr_exponent * pow(1.0-VolFr_TOTAL, VolFr_exponent-1.0);
              if ( Tau <= VolFr_min_vacant ) {
                Tau = 0.0;
                Tau__dcc =
                Tau__dfb = 0.0;
              }
            }

          const Real delta_RT = delta * (1.0 - exp(-delta_RT_a*RT-delta_RT_b*RT*RT));
          const Real lambda_RT = lambda * (RT/lambda_RT_r);

          for (std::size_t i=0; i<n_var_dofs; i++)
            {
              // RHS contribution
              Fe_var[1](i) += JxW[qp]*(
                                        cc_old * phi[i][qp] // capacity term
                                      + DT_2*( // source, sink terms
                                               kappa * Tau * (cc_old-cc_old*cc_old) * phi[i][qp]
                                             - delta_RT * cc_old * phi[i][qp]
                                             )
                                      );
              // RHS contribution
              Fe_var[2](i) += JxW[qp]*(
                                        fb_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               lambda_RT * Tau * fb_old * phi[i][qp]
                                             - diffusion * Tau * (GRAD_fb_old * dphi[i][qp])
                                             - haptotaxis * Tau * (GRAD_HU_old * fb_old * dphi[i][qp])
                                             )
                                      );

              for (std::size_t j=0; j<n_var_dofs; j++)
                {
                  // Matrix contribution
                  Ke_var[1][1](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        kappa * Tau__dcc * (cc_old-cc_old*cc_old) * phi[j][qp] * phi[i][qp]
                                                      + kappa * Tau * (1.0-2.0*cc_old) * phi[j][qp] * phi[i][qp]
                                                      - delta_RT * phi[j][qp] * phi[i][qp]
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
                                                      )
                                               );
                  Ke_var[1][2](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        kappa * Tau__dfb * (cc_old-cc_old*cc_old) * phi[j][qp] * phi[i][qp]
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[2][0](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                      - haptotaxis * Tau * (dphi[j][qp] * fb_old * dphi[i][qp])
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][1](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        lambda_RT * Tau__dcc * fb_old * phi[j][qp] * phi[i][qp]
                                                      - diffusion * Tau__dcc * phi[j][qp] * (GRAD_fb_old * dphi[i][qp])
                                                      - haptotaxis * Tau__dcc * phi[j][qp] * (GRAD_HU_old * fb_old * dphi[i][qp])
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][2](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        lambda_RT * Tau__dfb * fb_old * phi[j][qp] * phi[i][qp]
                                                      + lambda_RT * Tau * phi[j][qp] * phi[i][qp]
                                                      - diffusion * Tau__dfb * phi[j][qp] * (GRAD_fb_old * dphi[i][qp])
                                                      - diffusion * Tau * (dphi[j][qp] * dphi[i][qp])
                                                      - haptotaxis * Tau__dfb * phi[j][qp] * (GRAD_HU_old * fb_old * dphi[i][qp])
                                                      - haptotaxis * Tau * (GRAD_HU_old * phi[j][qp] * dphi[i][qp])
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
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
    es.get_system<TransientLinearImplicitSystem>("RIPF");
  libmesh_assert_equal_to(system.n_vars(), 3);

  ExplicitSystem & RT_system =
    es.get_system<ExplicitSystem>("RT");
  libmesh_assert_equal_to(RT_system.n_vars(), 3);

  std::vector<Number> soln;
  system.update_global_solution(soln);

  std::vector<Number> RT_soln;
  RT_system.update_global_solution(RT_soln);

  const Real HU_min = es.parameters.get<Real>("HU/min"),
             HU_max = es.parameters.get<Real>("HU/max");
  const Real RT_broad_frac = es.parameters.get<int>("RT_dose/broad/fractions"),
             RT_focus_frac = es.parameters.get<int>("RT_dose/focus/fractions"),
             RT_total_frac = RT_broad_frac + RT_focus_frac;
  // simulation time is expressed in "days"
  const int day = std::floor( system.time );

  for (const auto & node : mesh.node_ptr_range())
    {
      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );

      Real HU_, cc_, fb_;
      HU_ = soln[idof[0]]; if (HU_<HU_min) HU_ = HU_min; else if (HU_>HU_max) HU_ = HU_max;
      cc_ = soln[idof[1]]; if (cc_<0.0) cc_ = 0.0;
      fb_ = soln[idof[2]]; if (fb_<0.0) fb_ = 0.0;

      system.solution->set(idof[0], HU_);
      system.solution->set(idof[1], cc_);
      system.solution->set(idof[2], fb_);

      const dof_id_type RT_idof[] = { node->dof_number(RT_system.number(), 0, 0) ,
                                      node->dof_number(RT_system.number(), 1, 0) ,
                                      node->dof_number(RT_system.number(), 2, 0) };
      libmesh_assert( node->n_comp(RT_system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(RT_system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(RT_system.number(), 2) == 1 );

      const Real RT_broad_ = RT_soln[RT_idof[0]],
                 RT_focus_ = RT_soln[RT_idof[1]];
      Real RT_total_ = 0.0; // radiation therapy (RT) dose per fraction
      if      ( day < RT_broad_frac ) RT_total_ = RT_broad_ / RT_broad_frac * (day+1);
      else if ( day < RT_total_frac ) RT_total_ = RT_focus_ / RT_focus_frac * ((day+1)-RT_broad_frac) + RT_broad_;
      else                            RT_total_ = RT_broad_ + RT_focus_;

      RT_system.solution->set(RT_idof[2], RT_total_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  RT_system.solution->close();
  RT_system.update();
  // ...done
}
