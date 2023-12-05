#include "./solid_system.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

static void input (const std::string & , EquationSystems & );
static void initial_hcc (EquationSystems & , const std::string & );
static void assemble_hcc (EquationSystems & , const std::string & );
static void initial_fibres (EquationSystems & , const std::string & );
static void check_solution (EquationSystems & );
static void adaptive_remeshing (EquationSystems & , MeshRefinement &);

extern PerfLog plog;
static Parallel::Communicator * pm_ptr = 0;

void coupled_hcc (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);
  MeshRefinement amr(mesh);

  pm_ptr = & init.comm();

  // threaded assembly doesn't currently work with the moving mesh code
  if (libMesh::n_threads()>1) return;

  input("input.dat", es);

  // create the reaction-diffusion system
  TransientLinearImplicitSystem & model_rds =
    es.add_system<TransientLinearImplicitSystem>("HCC");
  model_rds.add_variable("l", FIRST, LAGRANGE);
  model_rds.add_variable("c", FIRST, LAGRANGE);
  model_rds.add_variable("n", FIRST, LAGRANGE);
  model_rds.attach_init_function(initial_hcc);
  model_rds.attach_assemble_function(assemble_hcc);

  // create the solid biomechanics system
  SolidSystem & model_sb = es.add_system<SolidSystem>("SolidSystem");
  model_sb.add_variable("x", FIRST, LAGRANGE);
  model_sb.add_variable("y", FIRST, LAGRANGE);
  model_sb.add_variable("z", FIRST, LAGRANGE);

  // create an auxiliary system that stores the initial mesh
  TransientExplicitSystem & aux_sys = es.add_system<TransientExplicitSystem>("SolidSystem::auxiliary");
  aux_sys.add_variable("undeformed_x", FIRST, LAGRANGE);
  aux_sys.add_variable("undeformed_y", FIRST, LAGRANGE);
  aux_sys.add_variable("undeformed_z", FIRST, LAGRANGE);

  // create an additional system for the displacement vector field
  ExplicitSystem & disp_sys = es.add_system<ExplicitSystem>("SolidSystem::displacement");
  disp_sys.add_variable("u_x", FIRST, LAGRANGE);
  disp_sys.add_variable("u_y", FIRST, LAGRANGE);
  disp_sys.add_variable("u_z", FIRST, LAGRANGE);

  ExplicitSystem & fibre_sys = es.add_system<ExplicitSystem>("SolidSystem::fibre");
  fibre_sys.add_variable("fibre_reference_x", CONSTANT, MONOMIAL);
  fibre_sys.add_variable("fibre_reference_y", CONSTANT, MONOMIAL);
  fibre_sys.add_variable("fibre_reference_z", CONSTANT, MONOMIAL);
  fibre_sys.add_variable("fibre_current_x", CONSTANT, MONOMIAL);
  fibre_sys.add_variable("fibre_current_y", CONSTANT, MONOMIAL);
  fibre_sys.add_variable("fibre_current_z", CONSTANT, MONOMIAL);
  fibre_sys.attach_init_function(initial_fibres);

  // create an additional system for the hydrostatic pressure (mean solid stress)
  ExplicitSystem & press_sys = es.add_system<ExplicitSystem>("SolidSystem::pressure");
  press_sys.add_variable("p", CONSTANT, MONOMIAL);

  // create an additional system for the Von Mises stress
  ExplicitSystem & von_mises_sys = es.add_system<ExplicitSystem>("SolidSystem::von_mises");
  von_mises_sys.add_variable("VM", CONSTANT, MONOMIAL);

  GmshIO(mesh).read(es.parameters.get<std::string>("input_GMSH"));
  mesh.prepare_for_use(es.parameters.get<bool>("mesh/skip_renumber_nodes_and_elements"));
  mesh.print_info();
  GmshIO(mesh).write(es.parameters.get<std::string>("output_GMSH"));
  es.init();
  es.print_info();

  // perform important initializations to the solver and associated systems
  model_sb.save_initial_mesh();

  Paraview_IO paraview(mesh);
  paraview.open_pvd(es.parameters.get<std::string>("output_PARAVIEW"));

  // save initial solution
  paraview.update_pvd(es);

  const std::set<int> ltp = export_integers(es.parameters.get<std::string>("loading_time_points"));
  const std::set<int> otp = export_integers(es.parameters.get<std::string>("output_time_points"));
  const std::set<int> rtp = export_integers(es.parameters.get<std::string>("remeshing_time_points"));

  es.parameters.set<Real> ("time") = 0.0;
  es.parameters.set<Real>("pseudo_time") = 0.0;

  const int n_time_step = es.parameters.get<int>("number_of_time_steps");
  for (int t=1; t<=n_time_step; t++)
    {
      // update the simulation time
      es.parameters.set<Real>("time") += es.parameters.get<Real>("time_step");
      model_rds.time = es.parameters.get<Real>("time");
      if (ltp.end()!=ltp.find(t))
        es.parameters.set<Real>("pseudo_time") += model_sb.deltat;

      libMesh::out << std::endl
                   << " ==== Step " << std::setw(4) << t << " out of " << std::setw(4) << n_time_step
                   << " (time=" << es.parameters.get<Real>("time") << ") ==== "
                   << std::endl;

      // update the solution (containers) for up to 2 steps behind
      *(model_rds.older_local_solution) = *(model_rds.old_local_solution);
      *(model_rds.old_local_solution) = *(model_rds.current_local_solution);
      // now solve the AD progression model
      model_rds.solve();

      // check the reaction-diffusion system solution is in order
      check_solution(es);

      // solve for the solid (mechanics) equilibrium
      if (ltp.end()!=ltp.find(t))
        model_sb.run_solver();

      // perform post-processing of the solid system
      if (ltp.end()!=ltp.find(t))
        model_sb.post_process();

      // advance the Newmark time solver and then
      // update the auxiliary system only
      if (ltp.end()!=ltp.find(t))
        model_sb.update_data();

      // check if to perform any adaptive mesh refinement/coarsening
      if (rtp.end()!=rtp.find(t))
        adaptive_remeshing(es, amr);

      // save current solution
      if (otp.end()!=otp.find(t))
        paraview.update_pvd(es, t);
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

  // remove if directory to store simulation results exists already
  if (0==global_processor_id())
    std::system(std::string("rm -rf "+DIR).c_str());
  pm_ptr->barrier();
  // create a directory to store simulation results
  if (0==global_processor_id())
    std::system(std::string("mkdir "+DIR).c_str());
  pm_ptr->barrier();
  // create a copy of the input file containing all model parameters
  if (0==global_processor_id())
    std::system(std::string("cp "+file_name+" "+DIR+file_name).c_str());
  pm_ptr->barrier();

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
  name = "output_PARAVIEW";
  es.parameters.set<std::string>(name) = DIR + in(name, "output4paraview");
  //
  name = "input_fibres";
  es.parameters.set<std::string>(name) = in(name, ".");
  if (0==global_processor_id() && name!=".")
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0);
  name = "number_of_time_steps";
  es.parameters.set<int>(name) = in(name, 1);
  name = "number_of_loading_steps";
  es.parameters.set<int>(name) = in(name, 1);
  name = "loading_step";
  es.parameters.set<Real>(name) =
    ( es.parameters.get<Real>("time_step") * es.parameters.get<int>("number_of_time_steps") )
    / es.parameters.get<int>("number_of_loading_steps");
  // verification checks
  libmesh_assert(es.parameters.get<int>("number_of_time_steps")>=
                 es.parameters.get<int>("number_of_loading_steps"));
  libmesh_assert(es.parameters.get<Real>("loading_step")>=
                 es.parameters.get<Real>("time_step"));

  if ( es.parameters.get<int>("number_of_time_steps") %
       es.parameters.get<int>("number_of_loading_steps") )
    {
      libmesh_error();
    }
  else
    {
      const int time2loading_step =
        es.parameters.get<int>("number_of_time_steps") /
        es.parameters.get<int>("number_of_loading_steps");
      //
      int t = time2loading_step;
      std::string s;
      while (t<=es.parameters.get<int>("number_of_time_steps"))
        {
          s += " " + std::to_string(t) + " ";
          t += time2loading_step;
        }
      //
      es.parameters.set<std::string>("loading_time_points") = s;
    }

  name = "output_step";
  es.parameters.set<int>(name) = in(name, 0);
  //
  if (0==es.parameters.get<int>("output_step"))
    {
      std::string s;
      s += " " + std::to_string(es.parameters.get<int>("number_of_time_steps")) + " ";
      //
      es.parameters.set<std::string>("output_time_points") = s;
    }
  else
    {
      int t = es.parameters.get<int>("output_step");
      std::string s;
      while (t<=es.parameters.get<int>("number_of_time_steps"))
        {
          s += " " + std::to_string(t) + " ";
          t += es.parameters.get<int>("output_step");
        }
      //
      es.parameters.set<std::string>("output_time_points") = s;
    }

  name = "remeshing_step";
  es.parameters.set<int>(name) = in(name, 0);
  //
  if (0==es.parameters.get<int>("remeshing_step"))
    {
      std::string s;
      s += " " + std::to_string(1+es.parameters.get<int>("number_of_time_steps")) + " ";
      //
      es.parameters.set<std::string>("remeshing_time_points") = s;
    }
  else
    {
      int t = es.parameters.get<int>("remeshing_step");
      std::string s;
      while (t<=es.parameters.get<int>("number_of_time_steps"))
        {
          s += " " + std::to_string(t) + " ";
          t += es.parameters.get<int>("remeshing_step");
        }
      //
      es.parameters.set<std::string>("remeshing_time_points") = s;
    }

  {
    name = "mesh/skip_renumber_nodes_and_elements";
    es.parameters.set<bool>(name) = in(name, true);
    //
    name = "mesh/AMR/max_steps";
    es.parameters.set<int>(name) = in(name, 0);
    name = "mesh/AMR/max_level";
    es.parameters.set<int>(name) = in(name, 3);
    name = "mesh/AMR/refine_percentage";
    es.parameters.set<Real>(name) = in(name, 0.5);
    name = "mesh/AMR/coarsen_percentage";
    es.parameters.set<Real>(name) = in(name, 0.5);
  }

  {
    name = "solver/quiet";
    es.parameters.set<bool>(name) = in(name, false);
    //
    name = "solver/nonlinear/max_nonlinear_iterations";
    es.parameters.set<int>(name) = in(name, 100);
    name = "solver/nonlinear/relative_step_tolerance";
    es.parameters.set<Real>(name) = in(name, 1.e-3);
    name = "solver/nonlinear/relative_residual_tolerance";
    es.parameters.set<Real>(name) = in(name, 1.e-8);
    name = "solver/nonlinear/absolute_residual_tolerance";
    es.parameters.set<Real>(name) = in(name, 1.e-8);
    name = "solver/nonlinear/require_reduction";
    es.parameters.set<bool>(name) = in(name, false);
    //
    name = "solver/linear/max_linear_iterations";
    es.parameters.set<int>(name) = in(name, 50000);
    name = "solver/linear/initial_linear_tolerance";
    es.parameters.set<Real>(name) = in(name, 1.e-3);
    //
    name = "solver/assembly_use_symmetry";
    es.parameters.set<bool>(name) = in(name, false);
  }

  name = "BCs";
  es.parameters.set<std::string>(name) = in(name, " 0 ");
  //
  const std::set<int> BCs = export_integers(es.parameters.get<std::string>(name));
  // read the Dirichlet boundary conditions
  for (auto bc : BCs)
    {
      Point displacement_vector;
      for (int d=0; d<3; d++)
        {
          name = "BC/"+std::to_string(bc)+"/displacement/"+std::to_string(d);
          displacement_vector(d) = in(name, 0.0);
        }
      std::cout << bc << ' ' << displacement_vector << std::endl;
      //
      name = "BC/"+std::to_string(bc)+"/displacement";
      es.parameters.set<Point>(name) = displacement_vector;
    }
  name = "BCs/displacement_penalty";
  es.parameters.set<Real>(name) = in(name, 1.e+5);

  name = "materials";
  es.parameters.set<std::string>(name) = in(name, " 0 ");
  //
  const std::set<int> materials = export_integers(es.parameters.get<std::string>(name));
  // read the materials (subdomain regions)
  for (auto m : materials)
    {
      name = "material/"+std::to_string(m)+"/Hyperelastic/Young";
      es.parameters.set<Real>(name) = in(name, 1.e+3);
      name = "material/"+std::to_string(m)+"/Hyperelastic/Poisson";
      es.parameters.set<Real>(name) = in(name, 0.3);
      name = "material/"+std::to_string(m)+"/Hyperelastic/FibreStiffness";
      es.parameters.set<Real>(name) = in(name, 0.0);
      name = "material/"+std::to_string(m)+"/Hyperelastic/VolumetricStretchRatio/rate_0";
      es.parameters.set<Real>(name) = in(name, 0.0);
      name = "material/"+std::to_string(m)+"/Hyperelastic/VolumetricStretchRatio/rate_1";
      es.parameters.set<Real>(name) = in(name, 0.0);
      name = "material/"+std::to_string(m)+"/Hyperelastic/VolumetricStretchRatio/rate_2";
      es.parameters.set<Real>(name) = in(name, 0.0);
    }

  {
    name = "range/active_tumor/min"; es.parameters.set<Real>(name) = in(name, 1.0e-12);
    name = "range/active_tumor/max"; es.parameters.set<Real>(name) = in(name, 1.0e+12);
    name = "range/necrotic/min"; es.parameters.set<Real>(name) = in(name, 1.0e-12);
    name = "range/necrotic/max"; es.parameters.set<Real>(name) = in(name, 1.0e+12);
    name = "range/total_cell/min"; es.parameters.set<Real>(name) = in(name, 1.0e-12);
    name = "range/total_cell/max"; es.parameters.set<Real>(name) = in(name, 1.0e+12);
    //
    name = "cells_min_capacity";
    es.parameters.set<Real>(name) = in(name, 0.0);
    name = "cells_max_capacity";
    es.parameters.set<Real>(name) = in(name, 1.0);
    name = "cells_max_capacity/exponent";
    es.parameters.set<Real>(name) = in(name, 1.0);
  }

  {
    // parameters for variable: l
    name = "produce/l"; es.parameters.set<Real>(name) = in(name, 0.);
    // parameters for variable: c
    name = "diffuse/c"; es.parameters.set<Real>(name) = in(name, 0.);
    name = "mechano/c"; es.parameters.set<Real>(name) = in(name, 0.);
    name = "produce/c"; es.parameters.set<Real>(name) = in(name, 0.);
    // parameters for variable: n
    name = "necrosis/l";        es.parameters.set<Real>(name) = in(name, 0.);
    name = "necrosis/c";        es.parameters.set<Real>(name) = in(name, 0.);
    name = "necrosis/pressure"; es.parameters.set<Real>(name) = in(name, 0.);
  }

  // ...done
}

void initial_hcc (EquationSystems & es,
                  const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "HCC");

  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  TransientLinearImplicitSystem & hcc_sys =
    es.get_system<TransientLinearImplicitSystem>("HCC");
  libmesh_assert_equal_to(hcc_sys.n_vars(), 3);

  hcc_sys.time = 0.0;

  std::ifstream fin(es.parameters.get<std::string>("input_nodal"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real l_, c_, n_;
      fin >> l_ >> c_ >> n_;

      const dof_id_type idof[] = { node->dof_number(hcc_sys.number(), 0, 0) ,
                                   node->dof_number(hcc_sys.number(), 1, 0) ,
                                   node->dof_number(hcc_sys.number(), 2, 0) };
      libmesh_assert( node->n_comp(hcc_sys.number(), 0) == 1 );
      libmesh_assert( node->n_comp(hcc_sys.number(), 1) == 1 );
      libmesh_assert( node->n_comp(hcc_sys.number(), 2) == 1 );

      hcc_sys.solution->set(idof[0], l_);
      hcc_sys.solution->set(idof[1], c_);
      hcc_sys.solution->set(idof[2], n_);
    }

  // close solution vector and update the system
  hcc_sys.solution->close();
  hcc_sys.update();
  // ...done
}

void assemble_hcc (EquationSystems & es,
                   const std::string & system_name)
{
  libmesh_ignore(es, system_name);
  libmesh_assert_equal_to(system_name, "HCC");

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & hcc_sys =
    es.get_system<TransientLinearImplicitSystem>("HCC");
  libmesh_assert_equal_to(hcc_sys.n_vars(), 3);

  FEType fe_type = hcc_sys.variable_type(0);

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
  const Real Kappa_k = es.parameters.get<Real>("cells_max_capacity");
  const Real ek = es.parameters.get<Real>("cells_max_capacity/exponent");
  const Real produce_l  = es.parameters.get<Real>("produce/l");
  const Real diffuse_c_ = es.parameters.get<Real>("diffuse/c"),
             mechano_c_ = es.parameters.get<Real>("mechano/c"),
             produce_c  = es.parameters.get<Real>("produce/c");
  const Real necrosis_l = es.parameters.get<Real>("necrosis/l") / Kappa_k,
             necrosis_c = es.parameters.get<Real>("necrosis/c") / Kappa_k,
             necrosis_P = es.parameters.get<Real>("necrosis/pressure") / Kappa_k;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      hcc_sys.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(3);
      for (unsigned int v=0; v<3; v++)
        hcc_sys.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

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

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Number l_old(0.0), c_old(0.0), n_old(0.0);
          Gradient GRAD_c_old({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              l_old += phi[l][qp] * hcc_sys.old_solution(dof_indices_var[0][l]);
              c_old += phi[l][qp] * hcc_sys.old_solution(dof_indices_var[1][l]);
              n_old += phi[l][qp] * hcc_sys.old_solution(dof_indices_var[2][l]);
              GRAD_c_old.add_scaled(dphi[l][qp], hcc_sys.old_solution(dof_indices_var[1][l]));
            }

          Gradient GRAD_sigma({0.0, 0.0, 0.0});

          Real Tau(0.0);
          Real Tau__dl(0.0), Tau__dc(0.0), Tau__dn(0.0);
          {
            const Real Te_ = (l_old+c_old+n_old) / Kappa_k;
            if (Te_<=0.0)
              {
                Tau = 1.0;
                Tau__dl = Tau__dc = Tau__dn =
                0.0;
              }
            else if (Te_>=1.0)
              {
                Tau = 0.0;
                Tau__dl = Tau__dc = Tau__dn =
                0.0;
              }
            else
              {
                Tau = pow(1.0-Te_, ek);
                Tau__dl = Tau__dc = Tau__dn =
                (-ek/Kappa_k) * pow(1.0-Te_, ek-1.0);
              }
          }

          const Real diffuse_c = (c_old>Lambda_k ? diffuse_c_ : 0.0),
                     mechano_c = (c_old>Lambda_k ? mechano_c_ : 0.0);

          for (std::size_t i=0; i<n_var_dofs; i++)
            {
              // RHS contribution
              Fe_var[0](i) += JxW[qp]*(
                                        l_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               produce_l * Tau * l_old * phi[i][qp]
                                             - necrosis_l * l_old * n_old * phi[i][qp]
                                             )
                                      );
              // RHS contribution
              Fe_var[1](i) += JxW[qp]*(
                                        c_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               produce_c * Tau * c_old * phi[i][qp]
                                             - necrosis_c * c_old * n_old * phi[i][qp]
                                             - diffuse_c * Tau * (GRAD_c_old * dphi[i][qp])
                                             - mechano_c * Tau * c_old * (GRAD_sigma * dphi[i][qp])
                                             )
                                      );
              // RHS contribution
              Fe_var[2](i) += JxW[qp]*(
                                        n_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               necrosis_l * l_old * n_old * phi[i][qp]
                                             + necrosis_c * c_old * n_old * phi[i][qp]
                                             )
                                      );

              for (std::size_t j=0; j<n_var_dofs; j++)
                {
                  // Matrix contribution
                  Ke_var[0][0](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        produce_l * Tau * phi[j][qp] * phi[i][qp]
                                                      + produce_l * Tau__dl * phi[j][qp] * l_old * phi[i][qp]
                                                      - necrosis_l * phi[j][qp] * n_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[0][1](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        produce_l * Tau__dc * phi[j][qp] * l_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[0][2](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        produce_l * Tau__dn * phi[j][qp] * l_old * phi[i][qp]
                                                      - necrosis_l * l_old * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[1][0](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        produce_c * Tau__dl * phi[j][qp] * c_old * phi[i][qp]
                                                      - diffuse_c * Tau__dl * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - mechano_c * Tau__dl * phi[j][qp] * c_old * (GRAD_sigma * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[1][1](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        produce_c * Tau * phi[j][qp] * phi[i][qp]
                                                      + produce_c * Tau__dc * phi[j][qp] * c_old * phi[i][qp]
                                                      - necrosis_c * phi[j][qp] * n_old * phi[i][qp]
                                                      - diffuse_c * Tau__dc * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - diffuse_c * Tau * (dphi[j][qp] * dphi[i][qp])
                                                      - mechano_c * Tau__dc * phi[j][qp] * c_old * (GRAD_sigma * dphi[i][qp])
                                                      - mechano_c * Tau * phi[j][qp] * (GRAD_sigma * dphi[i][qp])
                                                      )
                                               );
                  Ke_var[1][1](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        produce_c * Tau__dn * phi[j][qp] * c_old * phi[i][qp]
                                                      - necrosis_c * c_old * phi[j][qp] * phi[i][qp]
                                                      - diffuse_c * Tau__dn * phi[j][qp] * (GRAD_c_old * dphi[i][qp])
                                                      - mechano_c * Tau__dn * phi[j][qp] * c_old * (GRAD_sigma * dphi[i][qp])
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[2][0](i,j) += JxW[qp]*(
                                               - DT_2*(
                                                        necrosis_l * phi[j][qp] * n_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][1](i,j) += JxW[qp]*(
                                               - DT_2*(
                                                        necrosis_c * phi[j][qp] * n_old * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][2](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*(
                                                        necrosis_l * l_old * phi[j][qp] * phi[i][qp]
                                                      + necrosis_c * c_old * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                }
            }
        }

      hcc_sys.get_dof_map().constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
      //
      hcc_sys.get_system_matrix().add_matrix(Ke, dof_indices);
      hcc_sys.rhs->add_vector(Fe, dof_indices);
    }

  // ...done
}

void initial_fibres (EquationSystems & es,
                     const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "SolidSystem::fibre");

  const MeshBase& mesh = es.get_mesh();

  ExplicitSystem & fibre_sys = es.get_system<ExplicitSystem>("SolidSystem::fibre");
  libmesh_assert_equal_to(fibre_sys.n_vars(), 6);

  const std::string name(es.parameters.get<std::string>("input_fibres"));
  if ("."==name) return;

  std::ifstream fin(name.c_str());

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      Real x_, y_, z_;
      fin >> x_ >> y_ >> z_;
      const Real m = sqrt(pow2(x_)+pow2(y_)+pow2(z_));
      if (m<=1.0e-6) libmesh_error();

      std::vector<dof_id_type> dof_indices_F_var[6];
      for (unsigned int f=0; f<6; f++)
        fibre_sys.get_dof_map().dof_indices(elem, dof_indices_F_var[f], f);

      fibre_sys.current_local_solution->set(dof_indices_F_var[0][0], x_/m);
      fibre_sys.current_local_solution->set(dof_indices_F_var[1][0], y_/m);
      fibre_sys.current_local_solution->set(dof_indices_F_var[2][0], z_/m);

      fibre_sys.current_local_solution->set(dof_indices_F_var[3][0], x_/m);
      fibre_sys.current_local_solution->set(dof_indices_F_var[4][0], y_/m);
      fibre_sys.current_local_solution->set(dof_indices_F_var[5][0], z_/m);
    }

  // fill global solution vector from local ones
  fibre_sys.current_local_solution->close();
  (*fibre_sys.solution) = (*fibre_sys.current_local_solution);
  fibre_sys.solution->close();
  // update the system
  fibre_sys.update();
  // ...done
}

void check_solution (EquationSystems & es)
{
  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  TransientLinearImplicitSystem & hcc_sys =
    es.get_system<TransientLinearImplicitSystem>("HCC");
  libmesh_assert_equal_to(hcc_sys.n_vars(), 3);

  std::vector<Number> soln;
  hcc_sys.update_global_solution(soln);

  for (const auto & node : mesh.node_ptr_range())
    {
      const dof_id_type idof[] = { node->dof_number(hcc_sys.number(), 0, 0) ,
                                   node->dof_number(hcc_sys.number(), 1, 0) ,
                                   node->dof_number(hcc_sys.number(), 2, 0) };
      libmesh_assert( node->n_comp(hcc_sys.number(), 0) == 1 );
      libmesh_assert( node->n_comp(hcc_sys.number(), 1) == 1 );
      libmesh_assert( node->n_comp(hcc_sys.number(), 2) == 1 );

      Real l_, c_, n_;
      l_ = soln[idof[0]]; if (l_<0.0) l_ = 0.0;
      c_ = soln[idof[1]]; if (c_<0.0) c_ = 0.0;
      n_ = soln[idof[2]]; if (n_<0.0) n_ = 0.0;

      hcc_sys.solution->set(idof[0], l_);
      hcc_sys.solution->set(idof[1], c_);
      hcc_sys.solution->set(idof[2], n_);
    }

  // close solution vector and update the system
  hcc_sys.solution->close();
  hcc_sys.update();

  // ...done
}

void adaptive_remeshing (EquationSystems & es, MeshRefinement & amr)
{
  ExplicitSystem& press_sys =
    es.get_system<ExplicitSystem>("SolidSystem::pressure");
  libmesh_assert_equal_to(press_sys.n_vars(), 1);

  const unsigned int varno__p = press_sys.variable_number("p");

  ExplicitSystem & von_mises_sys =
    es.get_system<ExplicitSystem>("SolidSystem::von_mises");
  libmesh_assert_equal_to(von_mises_sys.n_vars(), 1);

  const unsigned int varno__VM = von_mises_sys.variable_number("VM");

  TransientLinearImplicitSystem & hcc_sys =
    es.get_system<TransientLinearImplicitSystem>("HCC");
  libmesh_assert_equal_to(hcc_sys.n_vars(), 3);

  const unsigned int varno__l = hcc_sys.variable_number("l"),
                     varno__c = hcc_sys.variable_number("c"),
                     varno__n = hcc_sys.variable_number("n");

  const Real refine_pct  = es.parameters.get<Real>("mesh/AMR/refine_percentage"),
             coarsen_pct = es.parameters.get<Real>("mesh/AMR/coarsen_percentage");

  const int max_steps = es.parameters.get<int>("mesh/AMR/max_steps"),
            max_level = es.parameters.get<int>("mesh/AMR/max_level");
  if (0 == max_steps) return;

  for (int r=0; r<max_steps; r++)
    {
      ErrorVector err_vector;
      ErrorEstimator::ErrorMap err_map;
      err_map.insert(std::make_pair(std::make_pair(&press_sys, varno__p), &err_vector));
      err_map.insert(std::make_pair(std::make_pair(&hcc_sys, varno__c), &err_vector));
      // evaluate the error for each active element using the estimator provided
      KellyErrorEstimator err_estimator;
      err_estimator.estimate_errors(es, err_map);
      // flag elements for coarsening / refinement by mean and std.
      amr.flag_elements_by_mean_stddev(err_vector, refine_pct, coarsen_pct, max_level);
      // refine / coarsen the flagged FEs
      amr.refine_and_coarsen_elements();
      // reinitializes the equations system object
      es.reinit();
    }

  // ...done
}
