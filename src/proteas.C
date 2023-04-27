#include "./utils.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

static void input (const std::string & , EquationSystems & );
static void initial_radiotherapy_dosage (EquationSystems & , const std::string & );
static void initial_hounsfield_unit (EquationSystems & , const std::string & );
static void initial_proteas_model (EquationSystems & , const std::string & );
static void check_solution (EquationSystems & );
static void adaptive_mesh_refinement (EquationSystems & , MeshRefinement & );

extern PerfLog plog;
static Parallel::Communicator * pm_ptr = 0;

void proteas (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);
  MeshRefinement amr(mesh);

  pm_ptr = & init.comm();

  input("input.dat", es);

  TransientLinearImplicitSystem & model =
    es.add_system<TransientLinearImplicitSystem>("PROTEAS_model");
  model.add_variable("nec", FIRST, LAGRANGE);
  model.add_variable("ter", FIRST, LAGRANGE);
  model.add_variable("oed", FIRST, LAGRANGE);
  model.add_variable("vsc", FIRST, LAGRANGE);
  model.add_variable("gmt", FIRST, LAGRANGE);
  model.add_variable("wmt", FIRST, LAGRANGE);
  model.attach_init_function(initial_proteas_model);

  ExplicitSystem & RTD =
    es.add_system<ExplicitSystem>("RTD");
  RTD.add_variable("RTD", FIRST, LAGRANGE);
  RTD.attach_init_function(initial_radiotherapy_dosage);

  ExplicitSystem & HU =
    es.add_system<ExplicitSystem>("HU");
  HU.add_variable("HU", FIRST, LAGRANGE);
  HU.attach_init_function(initial_hounsfield_unit);

  GmshIO(mesh).read(es.parameters.get<std::string>("input_GMSH"));
  mesh.prepare_for_use(es.parameters.get<bool>("mesh/skip_renumber_nodes_and_elements"));
  mesh.print_info();
  GmshIO(mesh).write(es.parameters.get<std::string>("output_GMSH"));
  es.init();
  es.print_info();

  Paraview_IO paraview(mesh);
  paraview.open_pvd(es.parameters.get<std::string>("output_PARAVIEW"));

  std::ofstream csv;
  if (0==global_processor_id())
    csv.open(es.parameters.get<std::string>("output_CSV"));

  // save initial solution
  paraview.update_pvd(es);

  const std::set<int> otp = export_integers(es.parameters.get<std::string>("output_time_points"));

  const int refinement_step = es.parameters.get<int>("refinement_step");
  const int n_t_step = es.parameters.get<int>("time_step_number");
  for (int t=1; t<=n_t_step; t++)
    {
      // update the simulation time
      es.parameters.set<Real>("time") += es.parameters.get<Real>("time_step");
      model.time = es.parameters.get<Real>("time");

      libMesh::out << " ==== Step " << std::setw(4) << t << " out of " << std::setw(4) << n_t_step
                   << " (Time=" << std::setw(9) << model.time << ") ==== "
                   << std::endl;

      // update the solution (containers) for up to 2 steps behind
      *(model.older_local_solution) = *(model.old_local_solution);
      *(model.old_local_solution) = *(model.current_local_solution);
      // now solve the AD progression model
      model.solve();

      check_solution(es);

      if (0 == t%refinement_step)
        adaptive_mesh_refinement(es, amr);

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
  name = "input_NodalData";
  es.parameters.set<std::string>(name) = in(name, "input.nd");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());
  //
  name = "input_RTDosage";
  es.parameters.set<std::string>(name) = in(name, "input.rtd");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());
  //
  name = "input_HU";
  es.parameters.set<std::string>(name) = in(name, "input.hu");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+es.parameters.get<std::string>(name)).c_str());
  //
  name = "output_Paraview";
  es.parameters.set<std::string>(name) = DIR + in(name, "output4paraview");
  //
  name = "output_CSV";
  es.parameters.set<std::string>(name) = DIR + in(name, "output.csv");

  es.parameters.set<Real>("time") = 0.0;

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = in(name, 1);
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

  // ...done
}

void initial_radiotherapy_dosage (EquationSystems & es,
                                  const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "RTD");

  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  ExplicitSystem & system =
    es.get_system<ExplicitSystem>("RTD");
  libmesh_assert_equal_to(system.n_vars(), 1);

  std::ifstream fin(es.parameters.get<std::string>("input_RTDosage"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real RTD_;
      fin >> RTD_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );

      system.solution->set(idof[0], RTD_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void initial_hounsfield_unit (EquationSystems & es,
                              const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "HU");

  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  ExplicitSystem & system =
    es.get_system<ExplicitSystem>("HU");
  libmesh_assert_equal_to(system.n_vars(), 1);

  std::ifstream fin(es.parameters.get<std::string>("input_HU"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real HU_;
      fin >> HU_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );

      system.solution->set(idof[0], HU_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void initial_proteas_model (EquationSystems & es,
                            const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "PROTEAS_model");

  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PROTEAS_model");
  libmesh_assert_equal_to(system.n_vars(), 6);

  es.parameters.set<Real> ("time") =
  system.time = 0.0;

  std::ifstream fin(es.parameters.get<std::string>("input_NodalData"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real nec_, ter_, oed_, vsc_, gmt_, wmt_;
      fin >> nec_ >> ter_ >> oed_ >> vsc_ >> gmt_ >> wmt_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) ,
                                   node->dof_number(system.number(), 3, 0) ,
                                   node->dof_number(system.number(), 4, 0) ,
                                   node->dof_number(system.number(), 5, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );
      libmesh_assert( node->n_comp(system.number(), 3) == 1 );
      libmesh_assert( node->n_comp(system.number(), 4) == 1 );
      libmesh_assert( node->n_comp(system.number(), 5) == 1 );

      system.solution->set(idof[0], nec_);
      system.solution->set(idof[1], ter_);
      system.solution->set(idof[2], oed_);
      system.solution->set(idof[3], vsc_);
      system.solution->set(idof[4], gmt_);
      system.solution->set(idof[5], wmt_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void check_solution (EquationSystems & es)
{
  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PROTEAS_model");
  libmesh_assert_equal_to(system.n_vars(), 6);

  std::vector<Number> soln;
  system.update_global_solution(soln);

  for (const auto & node : mesh.node_ptr_range())
    {
      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) ,
                                   node->dof_number(system.number(), 3, 0) ,
                                   node->dof_number(system.number(), 4, 0) ,
                                   node->dof_number(system.number(), 5, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );
      libmesh_assert( node->n_comp(system.number(), 3) == 1 );
      libmesh_assert( node->n_comp(system.number(), 4) == 1 );
      libmesh_assert( node->n_comp(system.number(), 5) == 1 );

      Real nec_, ter_, oed_, vsc_, gmt_, wmt_;
      nec_ = soln[idof[0]]; if (nec_<0.0) nec_ = 0.0;
      ter_ = soln[idof[1]]; if (ter_<0.0) ter_ = 0.0;
      oed_ = soln[idof[2]]; if (oed_<0.0) oed_ = 0.0;
      vsc_ = soln[idof[3]]; if (vsc_<0.0) vsc_ = 0.0;
      gmt_ = soln[idof[4]]; if (gmt_<0.0) gmt_ = 0.0;
      wmt_ = soln[idof[5]]; if (wmt_<0.0) wmt_ = 0.0;

      system.solution->set(idof[0], nec_);
      system.solution->set(idof[1], ter_);
      system.solution->set(idof[2], oed_);
      system.solution->set(idof[3], vsc_);
      system.solution->set(idof[4], gmt_);
      system.solution->set(idof[5], wmt_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void adaptive_mesh_refinement (EquationSystems & es, MeshRefinement & amr)
{
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PROTEAS_model");
  libmesh_assert_equal_to(system.n_vars(), 6);

  // ...done
}
