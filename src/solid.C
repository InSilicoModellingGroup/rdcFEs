#include "./solid_system.h"

static void input (const std::string & , EquationSystems & );

extern PerfLog plog;
static Parallel::Communicator * pm_ptr = 0;

void solid (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);

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
  aux_sys.add_variable("undefo_x", FIRST, LAGRANGE);
  aux_sys.add_variable("undefo_y", FIRST, LAGRANGE);
  aux_sys.add_variable("undefo_z", FIRST, LAGRANGE);

  // create an additional system for the displacement vector field
  ExplicitSystem & disp_sys = es.add_system<ExplicitSystem>("SolidSystem::displacement");
  disp_sys.add_variable("U_x", FIRST, LAGRANGE);
  disp_sys.add_variable("U_y", FIRST, LAGRANGE);
  disp_sys.add_variable("U_z", FIRST, LAGRANGE);

  // create an additional system for the hydrostatic pressure (mean solid stress)
  ExplicitSystem & press_sys = es.add_system<ExplicitSystem>("SolidSystem::pressure");
  press_sys.add_variable("p", CONSTANT, MONOMIAL);

  GmshIO(mesh).read(es.parameters.get<std::string>("input_GMSH"));
  mesh.prepare_for_use(es.parameters.get<bool>("mesh/skip_renumber_nodes_and_elements"));
  mesh.print_info();
  GmshIO(mesh).write(es.parameters.get<std::string>("output_GMSH"));
  es.init();
  es.print_info();

  // perform important initializations to the solver and associated systems
  model.save_initial_mesh();

  const std::string ex2_filename =
    es.parameters.get<std::string>("output_EXODUS");

  const unsigned int n_t_step = 1.0 / model.deltat + 1;
  for (unsigned int t_step=0; t_step<n_t_step; t_step++)
    {
      const Real progress = t_step * model.deltat;
      es.parameters.set<Real>("time") = progress;
      es.parameters.set<unsigned int>("step") = t_step;

      out << "\n==== Time Step " << std::setw(4) << t_step;
      out << " (" << std::fixed << std::setprecision(2) << std::setw(6) << (progress*100.) << "%)";
      out << ", time = " << std::setw(7) << model.time;
      out << " ====\n" << std::endl;

      // solve for the solid (mechanics) equilibrium
      model.solve();

      // fill global solution vector from local ones
      aux_sys.current_local_solution->close();
      (*aux_sys.solution) = (*aux_sys.current_local_solution);
      aux_sys.solution->close();

      // reinitialize all systems
      es.reinit();

      // perform post-processing of the solid system
      model.post_process();

      // advance the Newmark solver of the solid system
      model.time_solver->advance_timestep();

      // specifically update the auxiliary system only
      model.update_auxiliary();
    }

  ExodusII_IO(es.get_mesh()).write_equation_systems(ex2_filename, es);

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
  name = "output_EXODUS";
  es.parameters.set<std::string>(name) = DIR + in(name, "output.ex2");

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time";
  es.parameters.set<Real>(name) = 0.0;
  name = "phase";
  es.parameters.set<unsigned int>(name) = 0;

  name = "mesh/skip_renumber_nodes_and_elements";
  es.parameters.set<bool>(name) = in(name, true);

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
      name = "material/"+std::to_string(m)+"/Neohookean/Young";
      es.parameters.set<Real>(name) = in(name, 1.e+3);
      name = "material/"+std::to_string(m)+"/Neohookean/Poisson";
      es.parameters.set<Real>(name) = in(name, 0.3);
    }

  // ...done
}
