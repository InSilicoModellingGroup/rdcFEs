#include "./utils.h"

static void input (const std::string & , EquationSystems & );
static void initial_tracts (EquationSystems & , const std::string & );
static void initial_adpm (EquationSystems & , const std::string & );
static void assemble_adpm (EquationSystems & , const std::string & );
static void check_solution (EquationSystems & );

extern PerfLog plog;
static Parallel::Communicator * pm_ptr = 0;

void adpm (LibMeshInit & init)
{
  Mesh mesh(init.comm(), 3);
  EquationSystems es(mesh);

  pm_ptr = & init.comm();

  input("input.dat", es);

  TransientLinearImplicitSystem & model =
    es.add_system<TransientLinearImplicitSystem>("ADPM");
  model.add_variable("PrP", FIRST, LAGRANGE);
  model.add_variable("A_b", FIRST, LAGRANGE);
  model.add_variable("Tau", FIRST, LAGRANGE);
  model.attach_assemble_function(assemble_adpm);
  model.attach_init_function(initial_adpm);

  ExplicitSystem & tracts =
    es.add_system<ExplicitSystem>("Tracts");
  tracts.add_variable("TractX", CONSTANT, MONOMIAL);
  tracts.add_variable("TractY", CONSTANT, MONOMIAL);
  tracts.add_variable("TractZ", CONSTANT, MONOMIAL);
  tracts.attach_init_function(initial_tracts);

  GmshIO(mesh).read(es.parameters.get<std::string>("input_GMSH"));
  mesh.prepare_for_use(es.parameters.get<bool>("mesh/skip_renumber_nodes_and_elements"));
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

  es.parameters.set<RealVectorValue>("velocity") = RealVectorValue(0., 0., 0.);

  es.parameters.set<Real>("time") = 0.0;

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = in(name, 1);
  name = "output_step";
  es.parameters.set<int>(name) = in(name, 1);

  name = "mesh/skip_renumber_nodes_and_elements";
  es.parameters.set<bool>(name) = in(name, true);

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
      Real x_, y_, z_;
      fin >> x_ >> y_ >> z_;

      std::vector<std::vector<dof_id_type>> dof_indices_T_var(3);

      system.get_dof_map().dof_indices(elem, dof_indices_T_var[0], 0);
      system.solution->set(dof_indices_T_var[0][0], x_);
      system.get_dof_map().dof_indices(elem, dof_indices_T_var[1], 1);
      system.solution->set(dof_indices_T_var[1][0], y_);
      system.get_dof_map().dof_indices(elem, dof_indices_T_var[2], 2);
      system.solution->set(dof_indices_T_var[2][0], z_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
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
      Real PrP_, A_b_, Tau_;
      fin >> PrP_ >> A_b_ >> Tau_;

      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );

      system.solution->set(idof[0], PrP_);
      system.solution->set(idof[1], A_b_);
      system.solution->set(idof[2], Tau_);
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

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("ADPM");
  libmesh_assert_equal_to(system.n_vars(), 3);

  const System & tracts_system =
    es.get_system<System>("Tracts");
  libmesh_assert_equal_to(tracts_system.n_vars(), 3);

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
        tracts_system.get_dof_map().dof_indices(elem, dof_indices_T_var[l], l);

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
          Number PrP_old(0.0), A_b_old(0.0), Tau_old(0.0);
          Gradient GRAD_PrP_old({0.0, 0.0, 0.0}), GRAD_A_b_old({0.0, 0.0, 0.0}), GRAD_Tau_old({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              PrP_old += phi[l][qp] * system.old_solution(dof_indices_var[0][l]);
              A_b_old += phi[l][qp] * system.old_solution(dof_indices_var[1][l]);
              Tau_old += phi[l][qp] * system.old_solution(dof_indices_var[2][l]);
              GRAD_PrP_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[0][l]));
              GRAD_A_b_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[1][l]));
              GRAD_Tau_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[2][l]));
            }

          for (std::size_t i=0; i<n_var_dofs; i++)
            {
              // RHS contribution
              Fe_var[0](i) += JxW[qp]*(
                                        PrP_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                             - Pi_(A_b_old,transform_A_b) * PrP_old * phi[i][qp]
                                             - Pi_(Tau_old,transform_Tau) * PrP_old * phi[i][qp]
                                             - Pi_(PrP_old,decay_PrP) * PrP_old * phi[i][qp]
                                             //- (GRAD_PrP_old * velocity) * phi[i][qp]
                                             )
                                      );
              // RHS contribution
              Fe_var[1](i) += JxW[qp]*(
                                        A_b_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               SD_(A_b_old,produce_A_b) * A_b_old * phi[i][qp]
                                             + Pi_(A_b_old,transform_A_b) * PrP_old * phi[i][qp]
                                             - Pi_(A_b_old,decay_A_b) * A_b_old * phi[i][qp]
                                             - Pi_(A_b_old,diffuse_A_b) * (GRAD_A_b_old * dphi[i][qp])
                                             - Pi_(A_b_old,taxis_A_b) * (GRAD_A_b_old * tracts) * (tracts * dphi[i][qp])
                                             //- (GRAD_A_b_old * velocity) * phi[i][qp]
                                             )
                                      );
              // RHS contribution
              Fe_var[2](i) += JxW[qp]*(
                                        Tau_old * phi[i][qp] // capacity term
                                      + DT_2*( // transport, source, sink terms
                                               SD_(Tau_old,produce_Tau) * Tau_old * phi[i][qp]
                                             + Pi_(Tau_old,transform_Tau) * PrP_old * phi[i][qp]
                                             - Pi_(Tau_old,decay_Tau) * Tau_old * phi[i][qp]
                                             - Pi_(Tau_old,diffuse_Tau) * (GRAD_Tau_old * dphi[i][qp])
                                             - Pi_(Tau_old,taxis_Tau) * (GRAD_Tau_old * tracts) * (tracts * dphi[i][qp])
                                             //- (GRAD_Tau_old * velocity) * phi[i][qp]
                                             )
                                      );

              for (std::size_t j=0; j<n_var_dofs; j++)
                {
                  // Matrix contribution
                  Ke_var[0][0](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                      - Pi_(A_b_old,transform_A_b) * phi[j][qp] * phi[i][qp]
                                                      - Pi_(Tau_old,transform_Tau) * phi[j][qp] * phi[i][qp]
                                                      - Pi_(PrP_old,decay_PrP) * phi[j][qp] * phi[i][qp]
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[1][0](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        Pi_(A_b_old,transform_A_b) * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[1][1](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        SD_(A_b_old,produce_A_b) * phi[j][qp] * phi[i][qp]
                                                      - Pi_(A_b_old,decay_A_b) * phi[j][qp] * phi[i][qp]
                                                      - Pi_(A_b_old,diffuse_A_b) * (dphi[j][qp] * dphi[i][qp])
                                                      - Pi_(A_b_old,taxis_A_b) * (dphi[j][qp] * tracts) * (tracts * dphi[i][qp])
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
                                                      )
                                               );
                  // Matrix contribution
                  Ke_var[2][0](i,j) += JxW[qp]*(
                                               - DT_2*( // transport, source, sink terms
                                                        Pi_(Tau_old,transform_Tau) * phi[j][qp] * phi[i][qp]
                                                      )
                                               );
                  Ke_var[2][2](i,j) += JxW[qp]*(
                                                 phi[j][qp] * phi[i][qp] // capacity term
                                               - DT_2*( // transport, source, sink terms
                                                        SD_(Tau_old,produce_Tau) * phi[j][qp] * phi[i][qp]
                                                      - Pi_(Tau_old,decay_Tau) * phi[j][qp] * phi[i][qp]
                                                      - Pi_(Tau_old,diffuse_Tau) * (dphi[j][qp] * dphi[i][qp])
                                                      - Pi_(Tau_old,taxis_Tau) * (dphi[j][qp] * tracts) * (tracts * dphi[i][qp])
                                                      //- (dphi[j][qp] * velocity) * phi[i][qp]
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

                Number  A_b_old(0.0), Tau_old(0.0);
                Gradient GRAD_A_b_old({0.0, 0.0, 0.0}), GRAD_Tau_old({0.0, 0.0, 0.0});
                for (std::size_t l=0; l<n_var_dofs; l++)
                  {
                    A_b_old += psi[l][qp] * system.old_solution(dof_indices_var[1][l]);
                    Tau_old += psi[l][qp] * system.old_solution(dof_indices_var[2][l]);
                    GRAD_A_b_old.add_scaled(dpsi[l][qp], system.old_solution(dof_indices_var[1][l]));
                    GRAD_Tau_old.add_scaled(dpsi[l][qp], system.old_solution(dof_indices_var[2][l]));
                  }

                Number field[] = {0.0, 0.0, 0.0};

                field[1] = qface_normals[qp] * ( //velocity * A_b_old +
                                                 Pi_(A_b_old,diffuse_A_b) * GRAD_A_b_old );
                field[2] = qface_normals[qp] * ( //velocity * Tau_old +
                                                 Pi_(Tau_old,diffuse_Tau) * GRAD_Tau_old );

                for (std::size_t i=0; i<psi.size(); i++)
                  {
                    // RHS contribution
                    Fe_var[1](i) += JxW_face[qp] * (psi[i][qp]*field[1])
                                  * penalty;
                    Fe_var[2](i) += JxW_face[qp] * (psi[i][qp]*field[2])
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

void check_solution (EquationSystems & es)
{
  const MeshBase& mesh = es.get_mesh();

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("ADPM");
  libmesh_assert_equal_to(system.n_vars(), 3);

  std::vector<Number> soln;
  system.update_global_solution(soln);

  for (const auto & node : mesh.node_ptr_range())
    {
      const dof_id_type idof[] = { node->dof_number(system.number(), 0, 0) ,
                                   node->dof_number(system.number(), 1, 0) ,
                                   node->dof_number(system.number(), 2, 0) };
      libmesh_assert( node->n_comp(system.number(), 0) == 1 );
      libmesh_assert( node->n_comp(system.number(), 1) == 1 );
      libmesh_assert( node->n_comp(system.number(), 2) == 1 );

      Real PrP_, A_b_, Tau_;
      PrP_ = soln[idof[0]]; if (PrP_<0.0) PrP_ = 0.0;
      A_b_ = soln[idof[1]]; if (A_b_<0.0) A_b_ = 0.0;
      Tau_ = soln[idof[2]]; if (Tau_<0.0) Tau_ = 0.0;

      system.solution->set(idof[0], PrP_);
      system.solution->set(idof[1], A_b_);
      system.solution->set(idof[2], Tau_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}
