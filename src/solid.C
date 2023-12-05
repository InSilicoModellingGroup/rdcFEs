#include "./solid_system.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

static void input (const std::string & , EquationSystems & );
static void initial_fibres (EquationSystems & , const std::string & );
static void adaptive_mesh_refinement (EquationSystems & , MeshRefinement &);

extern PerfLog plog;
static Parallel::Communicator * pm_ptr = 0;

void solid (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);
  MeshRefinement amr(mesh);

  pm_ptr = & init.comm();

  // threaded assembly doesn't currently work with the moving mesh code
  if (libMesh::n_threads()>1) return;

  input("input.dat", es);

  // create the core (non-linear solid mechanics) system
  SolidSystem & model = es.add_system<SolidSystem>("SolidSystem");
  model.add_variable("x", FIRST, LAGRANGE);
  model.add_variable("y", FIRST, LAGRANGE);
  model.add_variable("z", FIRST, LAGRANGE);

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
  model.save_initial_mesh();

  Paraview_IO paraview(mesh);
  paraview.open_pvd(es.parameters.get<std::string>("output_PARAVIEW"));

  // save initial solution
  paraview.update_pvd(es);

  const std::set<int> otp = export_integers(es.parameters.get<std::string>("output_time_points"));

  const int refinement_step = es.parameters.get<int>("refinement_step");
  const int n_t_step = es.parameters.get<int>("time_step_number");
  for (int t=1; t<=n_t_step; t++)
    {
      const Real time = t * model.deltat;
      es.parameters.set<Real>("time") = time;
      es.parameters.set<unsigned int>("step") = t;

      libMesh::out << std::endl
                   << " ==== Step " << std::setw(4) << t << " out of " << std::setw(4) << n_t_step
                   << " (Time=" << std::setw(9) << time << ") ==== "
                   << std::endl;

      // solve for the solid (mechanics) equilibrium
      model.run_solver();

      // perform post-processing of the solid system
      model.post_process();

      // advance the Newmark solver of the solid system
      model.time_solver->advance_timestep();

      // specifically update the auxiliary system only
      model.update_auxiliary();

      if (0 == t%refinement_step)
        adaptive_mesh_refinement(es, amr);

      // save current solution
      if (otp.end()!=otp.find(t))
        {
          paraview.update_pvd(es, t);
        }
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
  name = "output_PARAVIEW";
  es.parameters.set<std::string>(name) = DIR + in(name, "output4paraview");
  //
  name = "input_fibres";
  es.parameters.set<std::string>(name) = in(name, ".");
  if (0==global_processor_id() && name!=".")
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = 1.0/es.parameters.get<Real>("time_step");

  name = "output_step";
  es.parameters.set<int>(name) = in(name, 0);
  name = "refinement_step";
  es.parameters.set<int>(name) = in(name, 1+es.parameters.get<int>("time_step_number"));

  std::string otp;
  if (0==es.parameters.get<int>("output_step"))
    {
      name = "output_time_points";
      otp = in(name, std::to_string(es.parameters.get<int>("time_step_number")));
      es.parameters.set<std::string>(name) = otp;
    }
  else
    {
      int t = es.parameters.get<int>("output_step");
      std::string otp;
      while (t<=es.parameters.get<int>("time_step_number"))
        {
          otp += " " + std::to_string(t) + " ";
          t += es.parameters.get<int>("output_step");
        }
      name = "output_time_points";
      es.parameters.set<std::string>(name) = otp;
    }

  name = "time";
  es.parameters.set<Real>(name) = 0.0;
  name = "phase";
  es.parameters.set<unsigned int>(name) = 0;

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

void adaptive_mesh_refinement (EquationSystems & es, MeshRefinement & amr)
{
  ExplicitSystem& press_sys =
    es.get_system<ExplicitSystem>("SolidSystem::pressure");
  libmesh_assert_equal_to(press_sys.n_vars(), 1);

  const unsigned int varno__p = press_sys.variable_number("p");

  ExplicitSystem & von_mises_sys =
    es.get_system<ExplicitSystem>("SolidSystem::von_mises");
  libmesh_assert_equal_to(von_mises_sys.n_vars(), 1);

  const unsigned int varno__VM = von_mises_sys.variable_number("VM");

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
      err_map.insert(std::make_pair(std::make_pair(&von_mises_sys, varno__VM), &err_vector));
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