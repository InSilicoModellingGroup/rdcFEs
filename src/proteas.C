#include "./utils.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

static void input (const std::string & , EquationSystems & );
static void initial_radiotherapy_dosage (EquationSystems & , const std::string & );
static void initial_hounsfield_unit (EquationSystems & , const std::string & );
static void initial_proteas_model (EquationSystems & , const std::string & );
static void assemble_proteas_model (EquationSystems & , const std::string & );
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
  model.add_variable("hos", FIRST, LAGRANGE); // host (healthy) cells
  model.add_variable("tum", FIRST, LAGRANGE); // tumour cells
  model.add_variable("nec", FIRST, LAGRANGE); // necrotic cells
  model.add_variable("vsc", FIRST, LAGRANGE); // vascular cells
  model.add_variable("oed", FIRST, LAGRANGE); // oedema
  model.attach_init_function(initial_proteas_model);
  model.attach_assemble_function(assemble_proteas_model);

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

  std::ifstream fin(es.parameters.get<std::string>("input_nodal_RTDosage"));

  const Real RT_sessions = es.parameters.get<int>("RT/sessions");
    
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

  std::ifstream fin(es.parameters.get<std::string>("input_nodal_HU"));

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
  libmesh_assert_equal_to(system.n_vars(), 5);

  es.parameters.set<Real> ("time") =
    system.time = 0.0;

  std::ifstream fin(es.parameters.get<std::string>("input_NodalData"));

  for (const auto & node : mesh.node_ptr_range())
    {
      Real hos_, tum_, nec_, vsc_, oed_;
      fin >> hos_ >> tum_ >> nec_ >> vsc_ >> oed_;

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

      system.solution->set(idof[0], hos_);
      system.solution->set(idof[1], tum_);
      system.solution->set(idof[2], nec_);
      system.solution->set(idof[3], vsc_);
      system.solution->set(idof[4], oed_);
    }

  // close solution vector and update the system
  system.solution->close();
  system.update();
  // ...done
}

void assemble_proteas_model (EquationSystems & es,
                             const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "PROTEAS_model");

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const int proteas_model_vars = 5;
  const int RTD_vars =1;
  const int HU_vars = 1;
  
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PROTEAS_model");
  libmesh_assert_equal_to(system.n_vars(), proteas_model_vars);

  const ExplicitSystem & RTD_system =
    es.get_system<ExplicitSystem>("RTD");
  libmesh_assert_equal_to(RTD_system.n_vars(), RTD_vars);

  const ExplicitSystem & HU_system =
    es.get_system<ExplicitSystem>("HU");
  libmesh_assert_equal_to(RTD_system.n_vars(), HU_vars);

  FEType fe_type = system.variable_type(0);
  FEType fe_type_RTD = RTD_system.variable_type(0);
  FEType fe_type_HU = HU_system.variable_type(0);

  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
  std::unique_ptr<FEBase> fe_RTD(FEBase::build(dim, fe_type_RTD));
  std::unique_ptr<FEBase> fe_HU(FEBase::build(dim, fe_type_HU));

  QGauss qrule(dim, fe_type.default_quadrature_order());

  fe->attach_quadrature_rule(&qrule);
  fe_RTD->attach_quadrature_rule(&qrule);
  fe_HU->attach_quadrature_rule(&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<Real>> & phi_RTD = fe_RTD->get_phi();
  const std::vector<std::vector<Real>> & phi_HU = fe_HU->get_phi();

  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  const std::vector<std::vector<RealGradient>> & dphi_RTD = fe_RTD->get_dphi();
  const std::vector<std::vector<RealGradient>> & dphi_HU = fe_HU->get_dphi();

  const Real DT_2 = es.parameters.get<Real>("time_step") / 2.0;

  const Real T_max = es.parameters.get<Real>("cells/total_capacity");

  const Real RT_max = es.parameters.get<Real>("radiotherapy/max_dosage");

  const Real rho_h = es.parameters.get<Real>("host/proliferation"),
    u_h = es.parameters.get<Real>("host/vsc_threshold"),
    delta_h = es.parameters.get<Real>("host/RT_death_rate"),
    a_RT_h = es.parameters.get<Real>("host/RT_exp_a"),
    b_RT_h = es.parameters.get<Real>("host/RT_exp_b"),
    alpha_n_h = es.parameters.get<Real>("host/necrosis_rate");

  const Real D_c = es.parameters.get<Real>("tumour/diffusion"),
    D_c_h = es.parameters.get<Real>("tumour/diffusion_host"),
    rho_c = es.parameters.get<Real>("tumour/proliferation"),
    u_c = es.parameters.get<Real>("tumour/vsc_threshold"),
    delta_c = es.parameters.get<Real>("tumour/RT_death_rate"),
    a_RT_c = es.parameters.get<Real>("tumour/RT_exp_a"),
    b_RT_c = es.parameters.get<Real>("tumour/RT_exp_b"),
    alpha_n_c = es.parameters.get<Real>("tumour/necrosis_rate");

  const Real iota_n = es.parameters.get<Real>("necrosis/clearance"),
    k_n = es.parameters.get<Real>("necrosis/slope"),
    u_n = es.parameters.get<Real>("necrosis/vsc_threshold");

  const Real rho_v = es.parameters.get<Real>("vascular/proliferation"),
    alpha_n_v = es.parameters.get<Real>("vascular/necrosis_rate");

  const Real D_e = es.parameters.get<Real>("oedema/diffusion"),
    rho_e = es.parameters.get<Real>("oedema/proliferation"),
    u_e = es.parameters.get<Real>("oedema/vsc_threshold"),
    e_e = es.parameters.get<Real>("oedema/oedema_threshold"),
    xi_e = es.parameters.get<Real>("oedema/RT_coeff"),
    p_RT_e = es.parameters.get<Real>("oedema/RT_exp"),
    a_e = es.parameters.get<Real>("oedema/reabsorption_rate");

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      system.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(proteas_model_vars);
      for (unsigned int v=0; v<proteas_model_vars; v++)
        system.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

      std::vector<std::vector<dof_id_type>> dof_indices_RTD_var(RTD_vars);
      for (unsigned int l=0; l<HU_vars; l++)
        RTD_system.get_dof_map().dof_indices(elem, dof_indices_RTD_var[l], l);

      std::vector<std::vector<dof_id_type>> dof_indices_HU_var(HU_vars);
      for (unsigned int l=0; l<HU_vars; l++)
        HU_system.get_dof_map().dof_indices(elem, dof_indices_HU_var[l], l);

      const unsigned int n_dofs     = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();
      const unsigned int n_RTD_var_dofs = dof_indices_RTD_var[0].size();
      const unsigned int n_HU_var_dofs = dof_indices_HU_var[0].size();

      DenseMatrix<Number> Ke(n_dofs, n_dofs);
      DenseSubMatrix<Number> Ke_var[proteas_model_vars][proteas_model_vars] =
	{
	  { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
	  { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
	  { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
	  { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) } ,
	  { DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke) }
	};
      
      for (unsigned int i=0; i<proteas_model_vars; i++)
        for (unsigned int j=0; j<proteas_model_vars; j++)
          Ke_var[i][j].reposition(i*n_var_dofs, j*n_var_dofs, n_var_dofs, n_var_dofs);

      DenseVector<Number> Fe(n_dofs);
      DenseSubVector<Number> Fe_var[proteas_model_vars] =
	{
	  DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)
	};
      for (unsigned int i=0; i<proteas_model_vars; i++)
        Fe_var[i].reposition(i*n_var_dofs, n_var_dofs);

      fe->reinit(elem);
      fe_RTD->reinit(elem);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          Number hos_old(0.0), tum_old(0.0), nec_old(0.0), vsc_old(0.0), oed_old(0.0);
          Gradient GRAD_hos_old({0.0, 0.0, 0.0}),GRAD_tum_old({0.0, 0.0, 0.0}),GRAD_oed_old({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_var_dofs; l++)
            {
              hos_old += phi[l][qp] * system.old_solution(dof_indices_var[0][l]);
              tum_old += phi[l][qp] * system.old_solution(dof_indices_var[1][l]);
              nec_old += phi[l][qp] * system.old_solution(dof_indices_var[2][l]);
	      vsc_old += phi[l][qp] * system.old_solution(dof_indices_var[3][l]);
	      oed_old += phi[l][qp] * system.old_solution(dof_indices_var[4][l]);
	      GRAD_hos_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[0][l]));
              GRAD_tum_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[1][l]));
	      GRAD_oed_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices_var[4][l]));
            }

          Number RTD(0.0);
          Gradient GRAD_RTD({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_RTD_var_dofs; l++)
            {
              RTD += phi_RTD[l][qp] * RTD_system.current_solution(dof_indices_RTD_var[0][l]);
              GRAD_RTD.add_scaled(dphi_RTD[l][qp], RTD_system.current_solution(dof_indices_RTD_var[0][l]));
            }
          // normalize the gradient vector
          {
            const Real l2norm = GRAD_RTD.norm();
            GRAD_RTD = l2norm ? GRAD_RTD.unit() : Gradient(0.0,0.0,0.0);
          }

          Number HU(0.0);
          Gradient GRAD_HU({0.0, 0.0, 0.0});
          for (std::size_t l=0; l<n_HU_var_dofs; l++)
            {
              HU += phi_HU[l][qp] * HU_system.current_solution(dof_indices_HU_var[0][l]);
              GRAD_HU.add_scaled(dphi_HU[l][qp], HU_system.current_solution(dof_indices_HU_var[0][l]));
            }
          // normalize the gradient vector
          {
            const Real l2norm = GRAD_HU.norm();
            GRAD_HU = l2norm ? GRAD_HU.unit() : Gradient(0.0,0.0,0.0);
          }

	  const Real T = hos_old + tum_old + nec_old + vsc_old;
	  const Real Kappa = 1 - T/T_max;
	  const Real dKappa = -1;
	  
	  const Real host_prol = rho_h * Kappa * heaviside(vsc_old - u_h);
	  const Real dhost_prol = rho_h * dKappa * heaviside(vsc_old - u_h);
	  const Real host_RT_death = delta_h * (1.0 - exp(- a_RT_h*RTD - b_RT_h*pow2(RTD)));
	  const Real host_nec = alpha_n_h * nec_old;

	  const Real tumour_prol = rho_c * Kappa * heaviside(vsc_old - u_c);
	  const Real dtumour_prol = rho_c * dKappa * heaviside(vsc_old - u_c);
	  const Real tumour_RT_death = delta_c * (1.0 - exp(- a_RT_c*RTD - b_RT_c*pow2(RTD)));
	  const Real tumour_nec = alpha_n_c * nec_old;
	  
	  const Real nec_prol = alpha_n_h * hos_old + alpha_n_c * tum_old + alpha_n_v * vsc_old;
	  const Real nec_clearance = iota_n*(1 - tanh(k_n*(vsc_old - u_n)));
	  const Real dnec_clearance_dv = iota_n* -k_n / (cosh(k_n*(vsc_old - u_n)) * cosh(k_n*(vsc_old - u_n)));

	  const Real vsc_prol = rho_v * Kappa * tum_old;
	  const Real dvsc_prol = rho_v * dKappa * tum_old;
	  const Real vsc_nec = alpha_n_v * nec_old;

	  const Real oed_prol = rho_e * tum_old * (1- tum_old);
	  const Real doed_prol_dc = rho_e * (1- 2*tum_old);
	  const Real oed_RT = xi_e * std::pow(RTD / RT_max,p_RT_e);
	  const Real oed_clearance = a_e * (1 - heaviside(vsc_old - u_e));

	  // source terms
          for (std::size_t i=0; i<n_var_dofs; i++)
            {
              // Host cells
              Fe_var[0](i) += JxW[qp]*(
				       hos_old * phi[i][qp] 
				       + DT_2*( 
					       + host_prol * hos_old * (1 - hos_old) * phi[i][qp]
					       - host_RT_death * hos_old * phi[i][qp]
					       - host_nec * hos_old * phi[i][qp]
						)
				       );
              // Tumour cells
              Fe_var[1](i) += JxW[qp]*(
				       tum_old * phi[i][qp] 
				       + DT_2*(
					       - D_c * Kappa * (GRAD_tum_old * dphi[i][qp])
					       - D_c_h * Kappa * (GRAD_hos_old * tum_old * dphi[i][qp])
					       + tumour_prol * tum_old * phi[i][qp]
					       - tumour_RT_death * tum_old * phi[i][qp]
					       - tumour_nec * tum_old * phi[i][qp]
					       )
				       );
              // Necrotic cells
              Fe_var[2](i) += JxW[qp]*(
				       nec_old * phi[i][qp] 
				       + DT_2*(
					       + nec_prol * nec_old * phi[i][qp]
					       - nec_clearance * nec_old * phi[i][qp]
					       )
				       );
	      // Vascular cells
              Fe_var[3](i) += JxW[qp]*(
				       vsc_old * phi[i][qp] 
				       + DT_2*( 
					       + vsc_prol * vsc_old * phi[i][qp]
					       - vsc_nec * vsc_old * phi[i][qp]
						)
				       );

	      // Oedema
              Fe_var[4](i) += JxW[qp]*(
				       oed_old * phi[i][qp] 
				       + DT_2*(
					       - D_e * (GRAD_oed_old * dphi[i][qp])
                                               + oed_prol * oed_old * phi[i][qp]
					       - oed_RT * oed_old * phi[i][qp]
					       - oed_clearance * oed_old * phi[i][qp]
					       )
				       );
	      
              for (std::size_t j=0; j<n_var_dofs; j++)
                {
                  // Matrix contribution
		  // hos
                  Ke_var[0][0](i,j) += JxW[qp]*(
						phi[j][qp] * phi[i][qp]
						- DT_2*(
							+ dhost_prol * hos_old * (1 - hos_old) * phi[j][qp] * phi[i][qp]
							+ host_prol * (1 - 2*hos_old) * phi[j][qp] *phi[i][qp]
							- host_RT_death * phi[j][qp] * phi[i][qp]
							- host_nec * phi[j][qp] * phi[i][qp]
							)
						);
                  Ke_var[0][1](i,j) += JxW[qp]*(
						- DT_2*( 
							+ dhost_prol * hos_old * (1 - hos_old) * phi[j][qp] * phi[i][qp] 
							 )
						);
                  Ke_var[0][2](i,j) += JxW[qp]*(
						- DT_2*(
							+ dhost_prol * hos_old * (1 - hos_old) * phi[j][qp] * phi[i][qp]
							- alpha_n_h * phi[j][qp] * hos_old * phi[i][qp]
							)
						);
                  Ke_var[0][3](i,j) += JxW[qp]*(
						- DT_2*(
							+ dhost_prol * hos_old * (1 - hos_old) * phi[j][qp] * phi[i][qp]
							)
						);
                  Ke_var[0][4](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
		  // tum
                  Ke_var[1][0](i,j) += JxW[qp]*(
						- DT_2*(
							- D_c * dKappa * phi[j][qp] * (GRAD_tum_old * dphi[i][qp])
							- D_c_h * dKappa * phi[j][qp] * (GRAD_hos_old * tum_old * dphi[i][qp])
							- D_c_h * Kappa * (dphi[j][qp] * tum_old * dphi[i][qp])
							+ dtumour_prol * phi[j][qp] * tum_old * phi[i][qp]
							)
						);
                  Ke_var[1][1](i,j) += JxW[qp]*(
						phi[j][qp] * phi[i][qp]
						- DT_2*(
							- D_c * dKappa * phi[j][qp] * (GRAD_tum_old * dphi[i][qp])
							- D_c * Kappa * (dphi[j][qp] * dphi[i][qp])
							+ dtumour_prol * phi[j][qp] * tum_old * phi[i][qp]
							+ tumour_prol * phi[j][qp] * phi[i][qp]
							- tumour_RT_death * phi[j][qp] * phi[i][qp]
							- tumour_nec * phi[j][qp] * phi[i][qp]
		
							)
						);
                  Ke_var[1][2](i,j) += JxW[qp]*(
						- DT_2*(
							- D_c * dKappa * phi[j][qp] * (GRAD_tum_old * dphi[i][qp])
							- D_c_h * dKappa * phi[j][qp] * (GRAD_hos_old * tum_old * dphi[i][qp])
							+ dtumour_prol * phi[j][qp] * tum_old * phi[i][qp]
							- alpha_n_c * phi[j][qp] * tum_old * phi[i][qp]
						       
							)
						);
                  Ke_var[1][3](i,j) += JxW[qp]*(
						- DT_2*(
							- D_c * dKappa * phi[j][qp] * (GRAD_tum_old * dphi[i][qp])
							- D_c_h * dKappa * phi[j][qp] * (GRAD_hos_old * tum_old * dphi[i][qp])
							+ dtumour_prol * phi[j][qp] * tum_old * phi[i][qp]
							)
						);
                  Ke_var[1][4](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
		  // nec
                  Ke_var[2][0](i,j) += JxW[qp]*(
						- DT_2*(
							+ alpha_n_h * phi[j][qp] * nec_old * phi[i][qp] )
						);
                  Ke_var[2][1](i,j) += JxW[qp]*(
						- DT_2*( 
							+ alpha_n_c * phi[j][qp] * nec_old * phi[i][qp] 
							 )
						);
                  Ke_var[2][2](i,j) += JxW[qp]*(
						phi[j][qp] * phi[i][qp]
						- DT_2*(
							+ nec_prol * phi[j][qp] * phi[i][qp]
							- nec_clearance * phi[j][qp] * phi[i][qp]
							)
						);
                  Ke_var[2][3](i,j) += JxW[qp]*(
						- DT_2*(
							+ alpha_n_v * phi[j][qp] * nec_old * phi[i][qp]
							- dnec_clearance_dv * phi[j][qp] * nec_old * phi[i][qp]
							)
						);
                  Ke_var[2][4](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
		  // vsc
		  Ke_var[3][0](i,j) += JxW[qp]*(
						- DT_2*(
							+ dvsc_prol * phi[j][qp] * vsc_old * phi[i][qp]
							)
						);
                  Ke_var[3][1](i,j) += JxW[qp]*(
						- DT_2*(
							+ dvsc_prol * phi[j][qp] * vsc_old * phi[i][qp]
							)
						);
                  Ke_var[3][2](i,j) += JxW[qp]*(
						- DT_2*(
							+ dvsc_prol * phi[j][qp] * vsc_old * phi[i][qp]
							- alpha_n_v * phi[j][qp] * vsc_old * phi[i][qp]
							)
						);
                  Ke_var[3][3](i,j) += JxW[qp]*(
						phi[j][qp] * phi[i][qp]
						- DT_2*(
							+ dvsc_prol * phi[j][qp] * vsc_old * phi[i][qp]
							+ vsc_prol * phi[j][qp] * phi[i][qp]  
							- vsc_nec * phi[j][qp] * phi[i][qp]
							)
						);
                  Ke_var[3][4](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
		  // oedema
		  Ke_var[4][0](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
                  Ke_var[4][1](i,j) += JxW[qp]*(
						- DT_2*( 
							+ doed_prol_dc * phi[j][qp] * oed_old * phi[i][qp] 
							 )
						);
                  Ke_var[4][2](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
                  Ke_var[4][3](i,j) += JxW[qp]*(
						- DT_2*(0)
						);
                  Ke_var[4][4](i,j) += JxW[qp]*(
						phi[j][qp] * phi[i][qp]
						- DT_2*(
							- D_e * (dphi[j][qp] * dphi[i][qp])
							+ oed_prol * phi[j][qp] * phi[i][qp]
							- oed_RT * phi[j][qp] * phi[i][qp]
							- oed_clearance * phi[j][qp] * phi[i][qp]
							)
						);
		}
	    }
	}
      system.get_dof_map().constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

      system.get_system_matrix().add_matrix(Ke, dof_indices);
      system.rhs->add_vector(Fe, dof_indices);
    }
}

void check_solution (EquationSystems & es)
{
  const MeshBase& mesh = es.get_mesh();
  libmesh_assert_equal_to(mesh.mesh_dimension(), 3);

  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("PROTEAS_model");
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

      Real hos_, tum_, nec_, vsc_, oed_;
      hos_ = soln[idof[0]]; if (hos_<0.0) hos_ = 0.0;
      tum_ = soln[idof[1]]; if (tum_<0.0) tum_ = 0.0;
      nec_ = soln[idof[2]]; if (nec_<0.0) nec_ = 0.0;
      vsc_ = soln[idof[3]]; if (vsc_<0.0) vsc_ = 0.0;
      oed_ = soln[idof[4]]; if (oed_<0.0) oed_ = 0.0;

      system.solution->set(idof[0], hos_);
      system.solution->set(idof[1], tum_);
      system.solution->set(idof[2], nec_);
      system.solution->set(idof[3], vsc_);
      system.solution->set(idof[4], oed_);
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
