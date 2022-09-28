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
