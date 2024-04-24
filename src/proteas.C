#include "./utils.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

static void input (const std::string & , EquationSystems & );
static void calc_capacity_matrix (EquationSystems & );
static void calc_rhs_vector (EquationSystems & );
static void explicit_solve (EquationSystems & );
static void initial_aux_data (EquationSystems & , const std::string & );
static void initial_proteas_model (EquationSystems & , const std::string & );
static void check_solution (EquationSystems & );
static void adaptive_remeshing (EquationSystems & , MeshRefinement & );

extern PerfLog plog;
static Parallel::Communicator * pm_ptr = 0;

void proteas (LibMeshInit & init, const std::string &inputFile)
{
  Mesh msh(init.comm(), 3);
  EquationSystems es(msh);
  MeshRefinement amr(msh);

  pm_ptr = & init.comm();

  input(inputFile, es);

  LinearImplicitSystem & sys_PROTEAS =
    es.add_system<LinearImplicitSystem>("PROTEAS");
  sys_PROTEAS.add_variable("hos", FIRST, LAGRANGE); // host (healthy) cells
  sys_PROTEAS.add_variable("tum", FIRST, LAGRANGE); // tumour cells
  sys_PROTEAS.add_variable("nec", FIRST, LAGRANGE); // necrotic cells
  sys_PROTEAS.add_variable("vsc", FIRST, LAGRANGE); // vascular cells
  sys_PROTEAS.add_variable("oed", FIRST, LAGRANGE); // oedema
  sys_PROTEAS.attach_init_function(initial_proteas_model);
  sys_PROTEAS.add_vector("rhs");
  sys_PROTEAS.add_vector("reciprocal_vector_of_diagonal_capacity_matrix");
  sys_PROTEAS.add_vector("previous_solution");

  ExplicitSystem & sys_PROTEAS_dt =
    es.add_system<ExplicitSystem>("PROTEAS (rate)");
  sys_PROTEAS_dt.add_variable("hos_dt", FIRST, LAGRANGE); // host (healthy) cells
  sys_PROTEAS_dt.add_variable("tum_dt", FIRST, LAGRANGE); // tumour cells
  sys_PROTEAS_dt.add_variable("nec_dt", FIRST, LAGRANGE); // necrotic cells
  sys_PROTEAS_dt.add_variable("vsc_dt", FIRST, LAGRANGE); // vascular cells
  sys_PROTEAS_dt.add_variable("oed_dt", FIRST, LAGRANGE); // oedema

  ExplicitSystem & sys_AUX =
    es.add_system<ExplicitSystem>("AUX");
  sys_AUX.add_variable("HU", FIRST, LAGRANGE);
  sys_AUX.add_variable("RTD", FIRST, LAGRANGE);
  sys_AUX.attach_init_function(initial_aux_data);

  GmshIO(msh).read(es.parameters.get<std::string>("input_GMSH"));
  msh.prepare_for_use(es.parameters.get<bool>("mesh/skip_renumber_nodes_and_elements"));
  msh.print_info();
  GmshIO(msh).write(es.parameters.get<std::string>("output_GMSH"));
  es.init();
  es.print_info();

  Paraview_IO paraview(msh);
  paraview.open_pvd(es.parameters.get<std::string>("output_Paraview"));

  std::ofstream csv;
  if (0==global_processor_id())
    csv.open(es.parameters.get<std::string>("output_CSV"));

  // save initial solution
  paraview.update_pvd(es);

  const std::set<int> otp = export_integers(es.parameters.get<std::string>("output_time_points"));

  Real & time = es.parameters.set<Real>("time");

  const int remeshing_step = es.parameters.get<int>("remeshing_step");
  const int n_t_step = es.parameters.get<int>("time_step_number");

  calc_capacity_matrix(es);

  for (int t=1; t<=n_t_step; t++)
    {
      // update the simulation time
      time += es.parameters.get<Real>("time_step");

      libMesh::out << " === Step " << std::setw(4) << t << "/" << std::setw(4) << n_t_step
                   << " (Time=" << std::setw(9) << time << ") === \r" << std::flush;

      calc_rhs_vector(es);

      explicit_solve(es);

      check_solution(es);

      if (0 == t%remeshing_step)
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
  const std::string DIR = in(name, DIR_default);

  // create a directory to store simulation results
  if (0==global_processor_id())
    std::system(std::string("mkdir -p "+DIR).c_str());
  // create a copy of the input file containing all model parameters
  if (0==global_processor_id())
    std::system(std::string("cp "+file_name+" "+DIR+"/input.dat").c_str());

  name = "input_GMSH";
  es.parameters.set<std::string>(name) = in(name, "input.msh");
  //
  name = "output_GMSH";
  es.parameters.set<std::string>(name) = (DIR+"/output.msh");
  //
  name = "input_nodal";
  es.parameters.set<std::string>(name) = in(name, "input.nd");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+"/input.nd").c_str());
  //
  name = "input_nodal_aux";
  es.parameters.set<std::string>(name) = in(name, "input_aux.nd");
  if (0==global_processor_id())
    std::system(std::string("cp "+es.parameters.get<std::string>(name)+" "+DIR+"/input_aux.nd").c_str());
  //
  name = "output_Paraview";
  es.parameters.set<std::string>(name) = (DIR+"/output.paraview");
  //
  name = "output_CSV";
  es.parameters.set<std::string>(name) = (DIR+"/output.csv");

  es.parameters.set<Real>("time") = 0.0;

  name = "time_step";
  es.parameters.set<Real>(name) = in(name, 1.0e-9);
  name = "time_step_number";
  es.parameters.set<int>(name) = in(name, 1);
  name = "output_step";
  es.parameters.set<int>(name) = in(name, 0);
  name = "remeshing_step";
  es.parameters.set<int>(name) = in(name, 1+es.parameters.get<int>("time_step_number"));

  name = "linear_solver_tolerance";
  es.parameters.set<Real>("linear solver tolerance") = in(name, 1.0e-6);
  name = "linear_solver_maximum_iterations";
  es.parameters.set<unsigned int>("linear solver maximum iterations") = in(name, 1000);

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

  {
    name = "cells/total_capacity"; es.parameters.set<Real>(name) = in(name, 1.0);

    name = "radiotherapy/min_dosage"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "radiotherapy/max_dosage"; es.parameters.set<Real>(name) = in(name, 1.0);

    name = "host/vsc_threshold"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "host/proliferation"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "host/RT_death_rate"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "host/RT_exp_a"; es.parameters.set<Real>(name) = in(name, 1.0);
    name = "host/RT_exp_b"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "host/necrosis_rate"; es.parameters.set<Real>(name) = in(name, 1.0);

    name = "tumour/vsc_threshold"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "tumour/proliferation"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "tumour/RT_death_rate"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "tumour/RT_exp_a"; es.parameters.set<Real>(name) = in(name, 1.0);
    name = "tumour/RT_exp_b"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "tumour/necrosis_rate"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "tumour/diffusion"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "tumour/diffusion_host"; es.parameters.set<Real>(name) = in(name, 0.0);

    name = "necrosis/vsc_threshold"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "necrosis/vsc_slope"; es.parameters.set<Real>(name) = in(name, 1.0);
    name = "necrosis/clearance_rate"; es.parameters.set<Real>(name) = in(name, 0.0);

    name = "vascular/proliferation"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "vascular/necrosis_rate"; es.parameters.set<Real>(name) = in(name, 0.0);

    name = "oedema/vsc_threshold"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "oedema/proliferation"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "oedema/RT_inflammation_rate"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "oedema/RT_exp"; es.parameters.set<Real>(name) = in(name, 1.0);
    name = "oedema/reabsorption_rate"; es.parameters.set<Real>(name) = in(name, 0.0);
    name = "oedema/diffusion"; es.parameters.set<Real>(name) = in(name, 0.0);
  }

  es.parameters.print();

  // ...done
}

void calc_capacity_matrix (EquationSystems & es)
{
  const MeshBase & mesh = es.get_mesh();
  libmesh_assert_equal_to(3, mesh.mesh_dimension());

  const unsigned int dim = mesh.mesh_dimension();
  const int MODEL_vars = 5;

  LinearImplicitSystem & sys_PROTEAS =
    es.get_system<LinearImplicitSystem>("PROTEAS");
  libmesh_assert_equal_to(sys_PROTEAS.n_vars(), MODEL_vars);

  FEType fe_type = sys_PROTEAS.variable_type(0);

  QGauss qrule(dim, fe_type.default_quadrature_order());

  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  const Real reciprocal_dt = 1.0 / es.parameters.get<Real>("time_step");

  // initialize the system matrix
  sys_PROTEAS.matrix->zero();

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      sys_PROTEAS.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(MODEL_vars);
      for (unsigned int v=0; v<MODEL_vars; v++)
        sys_PROTEAS.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

      const unsigned int n_dofs     = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      DenseMatrix<Number> Me(n_dofs, n_dofs);
      DenseSubMatrix<Number> Me_var[MODEL_vars][MODEL_vars] =
      {
        { DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me) } ,
        { DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me) } ,
        { DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me) } ,
        { DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me) } ,
        { DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me), DenseSubMatrix<Number>(Me) }
      };

      for (unsigned int i=0; i<MODEL_vars; i++)
        for (unsigned int j=0; j<MODEL_vars; j++)
          Me_var[i][j].reposition(i*n_var_dofs, j*n_var_dofs, n_var_dofs, n_var_dofs);

      fe->reinit(elem);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          for (unsigned int i=0; i<n_var_dofs; i++)
            {
              for (unsigned int j=0; j<n_var_dofs; j++)
                {
                  // Host (healthy) cells
                  Me_var[0][0](i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp])*reciprocal_dt;
                  // Tumour cells
                  Me_var[1][1](i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp])*reciprocal_dt;
                  // Necrotic cells
                  Me_var[2][2](i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp])*reciprocal_dt;
                  // Vascular cells
                  Me_var[3][3](i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp])*reciprocal_dt;
                  // Oedema
                  Me_var[4][4](i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp])*reciprocal_dt;
                }
            }
        }

      // row-sum method for lumping
      for (unsigned int a=0; a<n_dofs; a++)
        {
          Real Me_rowsum = 0.0;
          for (unsigned int b=0; b<n_dofs; b++)
            {
              Me_rowsum += Me(a,b);
              Me(a,b) = 0.0;
            }
          Me(a,a) = Me_rowsum;
        }

      sys_PROTEAS.get_dof_map().constrain_element_matrix(Me, dof_indices);

      sys_PROTEAS.matrix->add_matrix(Me, dof_indices);
    }

  NumericVector<Number> & diag_v =
    sys_PROTEAS.get_vector("reciprocal_vector_of_diagonal_capacity_matrix");
  diag_v.zero();

  // close the system matrix first
  sys_PROTEAS.matrix->close();
  // extract the diagonal of the lumped (diagonal) capacity matrix
  sys_PROTEAS.matrix->get_diagonal(diag_v);
  // evaluate the reciprocal vector
  diag_v.reciprocal();
  diag_v.close();

  // ...done
}

void calc_rhs_vector (EquationSystems & es)
{
  const MeshBase & mesh = es.get_mesh();
  libmesh_assert_equal_to(3, mesh.mesh_dimension());

  const unsigned int dim = mesh.mesh_dimension();
  const int MODEL_vars = 5;
  const int AUX_vars = 2;

  LinearImplicitSystem & sys_PROTEAS =
    es.get_system<LinearImplicitSystem>("PROTEAS");
  libmesh_assert_equal_to(sys_PROTEAS.n_vars(), MODEL_vars);

  const ExplicitSystem & sys_AUX =
    es.get_system<ExplicitSystem>("AUX");
  libmesh_assert_equal_to(sys_AUX.n_vars(), AUX_vars);

  FEType fe_type = sys_PROTEAS.variable_type(0);
  FEType fe_type_AUX = sys_AUX.variable_type(0);

  QGauss qrule(dim, fe_type.default_quadrature_order());

  std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
  fe->attach_quadrature_rule(&qrule);

  std::unique_ptr<FEBase> fe_AUX(FEBase::build(dim, fe_type_AUX));
  fe_AUX->attach_quadrature_rule(&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<Real>> & phi_AUX = fe_AUX->get_phi();

  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
  const std::vector<std::vector<RealGradient>> & dphi_AUX = fe_AUX->get_dphi();

  const Real T_max = es.parameters.get<Real>("cells/total_capacity");
  const Real RT_min = es.parameters.get<Real>("radiotherapy/min_dosage"),
             RT_max = es.parameters.get<Real>("radiotherapy/max_dosage");

  const Real vsc_h   = es.parameters.get<Real>("host/vsc_threshold"),
             rho_h   = es.parameters.get<Real>("host/proliferation"),
             delta_h = es.parameters.get<Real>("host/RT_death_rate"),
             a_RT_h  = es.parameters.get<Real>("host/RT_exp_a"),
             b_RT_h  = es.parameters.get<Real>("host/RT_exp_b"),
             nu_h    = es.parameters.get<Real>("host/necrosis_rate");

  const Real vsc_c   = es.parameters.get<Real>("tumour/vsc_threshold"),
             rho_c   = es.parameters.get<Real>("tumour/proliferation"),
             delta_c = es.parameters.get<Real>("tumour/RT_death_rate"),
             a_RT_c  = es.parameters.get<Real>("tumour/RT_exp_a"),
             b_RT_c  = es.parameters.get<Real>("tumour/RT_exp_b"),
             nu_c    = es.parameters.get<Real>("tumour/necrosis_rate");
  const Real D_c     = es.parameters.get<Real>("tumour/diffusion"),
             D_c_h   = es.parameters.get<Real>("tumour/diffusion_host");

  const Real vsc_n  = es.parameters.get<Real>("necrosis/vsc_threshold"),
             vsc_k  = es.parameters.get<Real>("necrosis/vsc_slope"),
             iota_n = es.parameters.get<Real>("necrosis/clearance_rate");

  const Real rho_v = es.parameters.get<Real>("vascular/proliferation"),
             nu_v  = es.parameters.get<Real>("vascular/necrosis_rate");

  const Real vsc_e  = es.parameters.get<Real>("oedema/vsc_threshold"),
             rho_e  = es.parameters.get<Real>("oedema/proliferation"),
             xi_e   = es.parameters.get<Real>("oedema/RT_inflammation_rate"),
             c_RT_e = es.parameters.get<Real>("oedema/RT_exp"),
             psi_e  = es.parameters.get<Real>("oedema/reabsorption_rate");
  const Real D_e    = es.parameters.get<Real>("oedema/diffusion");

  // initialize the system rhs vector
  sys_PROTEAS.rhs->zero();

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      sys_PROTEAS.get_dof_map().dof_indices(elem, dof_indices);

      std::vector<std::vector<dof_id_type>> dof_indices_var(MODEL_vars);
      for (unsigned int v=0; v<MODEL_vars; v++)
        sys_PROTEAS.get_dof_map().dof_indices(elem, dof_indices_var[v], v);

      std::vector<std::vector<dof_id_type>> dof_indices_AUX_var(AUX_vars);
      for (unsigned int l=0; l<AUX_vars; l++)
        sys_AUX.get_dof_map().dof_indices(elem, dof_indices_AUX_var[l], l);

      const unsigned int n_dofs     = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();
      const unsigned int n_AUX_var_dofs = dof_indices_AUX_var[0].size();

      DenseVector<Number> Fe(n_dofs);
      DenseSubVector<Number> Fe_var[MODEL_vars] =
      {
        DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)
      };
      for (unsigned int i=0; i<MODEL_vars; i++)
        Fe_var[i].reposition(i*n_var_dofs, n_var_dofs);

      fe->reinit(elem);
      fe_AUX->reinit(elem);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          Number hos(0.0), tum(0.0), nec(0.0), vsc(0.0), oed(0.0);
          Gradient GRAD_hos({0.0, 0.0, 0.0}), GRAD_tum({0.0, 0.0, 0.0}), GRAD_oed({0.0, 0.0, 0.0});
          for (unsigned int l=0; l<n_var_dofs; l++)
            {
              hos += phi[l][qp] * sys_PROTEAS.current_solution(dof_indices_var[0][l]);
              tum += phi[l][qp] * sys_PROTEAS.current_solution(dof_indices_var[1][l]);
              nec += phi[l][qp] * sys_PROTEAS.current_solution(dof_indices_var[2][l]);
              vsc += phi[l][qp] * sys_PROTEAS.current_solution(dof_indices_var[3][l]);
              oed += phi[l][qp] * sys_PROTEAS.current_solution(dof_indices_var[4][l]);
              GRAD_hos.add_scaled(dphi[l][qp], sys_PROTEAS.current_solution(dof_indices_var[0][l]));
              GRAD_tum.add_scaled(dphi[l][qp], sys_PROTEAS.current_solution(dof_indices_var[1][l]));
              GRAD_oed.add_scaled(dphi[l][qp], sys_PROTEAS.current_solution(dof_indices_var[4][l]));
            }
          {
            GRAD_hos = GRAD_hos.norm()>1.0e-6 ? GRAD_hos.unit() : Gradient(0.0, 0.0, 0.0);
          }

          Number HU(0.0), RTD(0.0);
          Gradient GRAD_HU({0.0, 0.0, 0.0}), GRAD_RTD({0.0, 0.0, 0.0});
          for (unsigned int l=0; l<n_AUX_var_dofs; l++)
            {
              HU  += phi_AUX[l][qp] * sys_AUX.current_solution(dof_indices_AUX_var[0][l]);
              RTD += phi_AUX[l][qp] * sys_AUX.current_solution(dof_indices_AUX_var[1][l]);
              GRAD_HU.add_scaled( dphi_AUX[l][qp], sys_AUX.current_solution(dof_indices_AUX_var[0][l]));
              GRAD_RTD.add_scaled(dphi_AUX[l][qp], sys_AUX.current_solution(dof_indices_AUX_var[1][l]));
            }
          {
            GRAD_HU  = GRAD_HU.norm() >1.0e-6 ? GRAD_HU.unit()  : Gradient(0.0, 0.0, 0.0);
            GRAD_RTD = GRAD_RTD.norm()>1.0e-6 ? GRAD_RTD.unit() : Gradient(0.0, 0.0, 0.0);
          }

          Real Kappa;
          {
            const Real T = hos + tum + nec + vsc;
            Kappa = 1.0 - (T/T_max);
            Kappa = std::min(std::max(Kappa,0.0),1.0);
          }

          Real Radio;
          if      (RTD<RT_min) Radio = 0.0;
          else if (RTD>RT_max) Radio = 1.0;
          else Radio = 1.0 - exp(-RTD*(a_RT_h+b_RT_h*RTD));

          const Real Omicron = std::pow(RTD/RT_max, c_RT_e);

          for (unsigned int i=0; i<n_var_dofs; i++)
            {
              // Host (healthy) cells
              Fe_var[0](i) += JxW[qp]*(
                                      //
                                        rho_h * heaviside(vsc-vsc_h) * Kappa * hos * (1.0-hos) * phi[i][qp]
                                      //
                                      - delta_h * Radio * hos * phi[i][qp]
                                      //
                                      - nu_h * nec * hos * phi[i][qp]
                                      );
              // Tumour cells
              Fe_var[1](i) += JxW[qp]*(
                                      //
                                        rho_c * heaviside(vsc-vsc_c) * Kappa * tum * phi[i][qp]
                                      //
                                      - delta_c * Radio * tum * phi[i][qp]
                                      //
                                      - nu_c * nec * tum * phi[i][qp]
                                      //
                                      - D_c * Kappa * (GRAD_tum * dphi[i][qp])
                                      //
                                      - D_c_h * Kappa * (GRAD_hos * tum * dphi[i][qp])
                                      );
              // Necrotic cells
              Fe_var[2](i) += JxW[qp]*(
                                      //
                                        nu_h * nec * hos * phi[i][qp]
                                      + nu_c * nec * tum * phi[i][qp]
                                      + nu_v * nec * vsc * phi[i][qp]
                                      //
                                      - iota_n * (1.0-tanh(vsc_k*(vsc-vsc_n))) * nec * phi[i][qp]
                                      );
              // Vascular cells
              Fe_var[3](i) += JxW[qp]*(
                                      //
                                        rho_v * Kappa * tum * vsc * phi[i][qp]
                                      //
                                      - nu_v * nec * vsc * phi[i][qp]
                                      );
              // Oedema
              Fe_var[4](i) += JxW[qp]*(
                                      //
                                        rho_e * tum * (1.0-tum) * oed * phi[i][qp]
                                      //
                                      - xi_e * Omicron * oed * phi[i][qp]
                                      //
                                      - psi_e * (1.0-heaviside(vsc-vsc_e)) * oed * phi[i][qp]
                                      //
                                      - D_e * (GRAD_oed * dphi[i][qp])
                                      );
            }
        }

      sys_PROTEAS.get_dof_map().constrain_element_vector(Fe, dof_indices);

      sys_PROTEAS.rhs->add_vector(Fe, dof_indices);
    }

  const NumericVector<Number> & diag_v =
    sys_PROTEAS.get_vector("reciprocal_vector_of_diagonal_capacity_matrix");

  // close the system rhs vector first
  sys_PROTEAS.rhs->close();
  // calculate the effective rhs vector
  sys_PROTEAS.get_vector("rhs").zero();
  sys_PROTEAS.get_vector("rhs").pointwise_mult(*sys_PROTEAS.rhs, diag_v);
  sys_PROTEAS.get_vector("rhs").close();

  // ...done
}

void explicit_solve (EquationSystems & es)
{
  const MeshBase & mesh = es.get_mesh();
  libmesh_assert_equal_to(3, mesh.mesh_dimension());

  const unsigned int dim = mesh.mesh_dimension();
  const int MODEL_vars = 5;

  LinearImplicitSystem & sys_PROTEAS =
    es.get_system<LinearImplicitSystem>("PROTEAS");
  libmesh_assert_equal_to(sys_PROTEAS.n_vars(), MODEL_vars);

  ExplicitSystem & sys_PROTEAS_dt =
    es.get_system<ExplicitSystem>("PROTEAS (rate)");

  const Real reciprocal_dt = 1.0 / es.parameters.get<Real>("time_step");

  // update the previous solution
  sys_PROTEAS.get_vector("previous_solution").zero();
  sys_PROTEAS.get_vector("previous_solution").add(*sys_PROTEAS.solution);
  sys_PROTEAS.get_vector("previous_solution").close();
  // update the current solution
  sys_PROTEAS.solution->add(sys_PROTEAS.get_vector("rhs"));
  sys_PROTEAS.solution->close();
  // update the entire system
  sys_PROTEAS.update();

  // update the rate solution
  sys_PROTEAS_dt.solution->zero();
  sys_PROTEAS_dt.solution->add( reciprocal_dt, *sys_PROTEAS.solution);
  sys_PROTEAS_dt.solution->add(-reciprocal_dt, sys_PROTEAS.get_vector("previous_solution"));
  sys_PROTEAS_dt.solution->close();
  // update the entire system
  sys_PROTEAS_dt.update();

  // ...done
}

void initial_aux_data (EquationSystems & es,
                       const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "AUX");

  const MeshBase & mesh = es.get_mesh();
  libmesh_assert_equal_to(3, mesh.mesh_dimension());

  const unsigned int dim = mesh.mesh_dimension();
  const int AUX_vars = 2;

  ExplicitSystem & sys_AUX =
    es.get_system<ExplicitSystem>("AUX");
  libmesh_assert_equal_to(sys_AUX.n_vars(), AUX_vars);

  std::ifstream fin(es.parameters.get<std::string>("input_nodal_aux"));
  if (! fin.is_open())
    {
      std::cout << "ERROR: Failed to open nodal aux input file "
                << es.parameters.get<std::string>("input_nodal_aux") << std::endl;
      exit(1);
    }

  for (const auto & node : mesh.node_ptr_range())
    {
      Real HU_, RTD_;

      std::string line;
      while ( std::getline(fin,line) )
        {
          // ignore empty lines and lines starting with '#'
          if (line.empty() || line[0] == '#') continue;
          // read all 2 species in consequtive order
          std::istringstream iss(line);
          if (iss >> HU_ >> RTD_)
            {
              break;
            }
          else
            {
              std::cout << "ERROR: Nodal input aux file failed to read line: " << line << std::endl;
              exit(1);
            }
          //
        }

      const dof_id_type idof[] = { node->dof_number(sys_AUX.number(), 0, 0),
                                   node->dof_number(sys_AUX.number(), 1, 0) };
      libmesh_assert( node->n_comp(sys_AUX.number(), 0) == 1 );
      libmesh_assert( node->n_comp(sys_AUX.number(), 1) == 1 );

      sys_AUX.solution->set(idof[0], HU_);
      sys_AUX.solution->set(idof[1], RTD_);
    }

  // close solution vector and update the system
  sys_AUX.solution->close();
  sys_AUX.update();

  // ...done
}

void initial_proteas_model (EquationSystems & es,
                            const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to(system_name, "PROTEAS");

  const MeshBase & mesh = es.get_mesh();
  libmesh_assert_equal_to(3, mesh.mesh_dimension());

  const unsigned int dim = mesh.mesh_dimension();
  const int MODEL_vars = 5;

  ExplicitSystem & sys_PROTEAS =
    es.get_system<ExplicitSystem>("PROTEAS");
  libmesh_assert_equal_to(sys_PROTEAS.n_vars(), MODEL_vars);

  std::ifstream fin(es.parameters.get<std::string>("input_nodal"));
  if (! fin.is_open())
    {
      std::cout << "ERROR: Failed to open nodal input file "
                << es.parameters.get<std::string>("input_nodal") << std::endl;
      exit(1);
    }

  for (const auto & node : mesh.node_ptr_range())
    {
      Real hos_, tum_, nec_, vsc_, oed_;

      std::string line;
      while (std::getline(fin,line))
        {
          // ignore empty lines and lines starting with '#'
          if (line.empty() || line[0] == '#') continue;
          // read all 5 species in consequtive order
          std::istringstream iss(line);
          if (iss >> hos_ >> tum_ >> nec_ >> vsc_ >> oed_)
            {
              break;
            }
          else
            {
              std::cout << "ERROR: Nodal input file failed to read line: " << line << std::endl;
              exit(1);
            }
          //
        }

      const dof_id_type idof[] = { node->dof_number(sys_PROTEAS.number(), 0, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 1, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 2, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 3, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 4, 0) };

      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 0) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 1) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 2) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 3) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 4) == 1 );

      sys_PROTEAS.solution->set(idof[0], hos_);
      sys_PROTEAS.solution->set(idof[1], tum_);
      sys_PROTEAS.solution->set(idof[2], nec_);
      sys_PROTEAS.solution->set(idof[3], vsc_);
      sys_PROTEAS.solution->set(idof[4], oed_);
    }

  fin.close();
  // close solution vector and update the system
  sys_PROTEAS.solution->close();
  sys_PROTEAS.update();

  // ...done
}

void check_solution (EquationSystems & es)
{
  const MeshBase & mesh = es.get_mesh();
  libmesh_assert_equal_to(3, mesh.mesh_dimension());

  const unsigned int dim = mesh.mesh_dimension();
  const int MODEL_vars = 5;

  ExplicitSystem & sys_PROTEAS =
    es.get_system<ExplicitSystem>("PROTEAS");
  libmesh_assert_equal_to(sys_PROTEAS.n_vars(), MODEL_vars);

  std::vector<Number> soln;
  sys_PROTEAS.update_global_solution(soln);

  for (const auto & node : mesh.node_ptr_range())
    {
      const dof_id_type idof[] = { node->dof_number(sys_PROTEAS.number(), 0, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 1, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 2, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 3, 0) ,
                                   node->dof_number(sys_PROTEAS.number(), 4, 0) };
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 0) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 1) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 2) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 3) == 1 );
      libmesh_assert( node->n_comp(sys_PROTEAS.number(), 4) == 1 );

      Real hos_, tum_, nec_, vsc_, oed_;
      hos_ = soln[idof[0]]; if (hos_<0.0) hos_ = 0.0;
      tum_ = soln[idof[1]]; if (tum_<0.0) tum_ = 0.0;
      nec_ = soln[idof[2]]; if (nec_<0.0) nec_ = 0.0;
      vsc_ = soln[idof[3]]; if (vsc_<0.0) vsc_ = 0.0;
      oed_ = soln[idof[4]]; if (oed_<0.0) oed_ = 0.0;

      sys_PROTEAS.solution->set(idof[0], hos_);
      sys_PROTEAS.solution->set(idof[1], tum_);
      sys_PROTEAS.solution->set(idof[2], nec_);
      sys_PROTEAS.solution->set(idof[3], vsc_);
      sys_PROTEAS.solution->set(idof[4], oed_);
    }

  // close solution vector and update the system
  sys_PROTEAS.solution->close();
  sys_PROTEAS.update();

  // ...done
}

void adaptive_remeshing (EquationSystems & es, MeshRefinement & amr)
{

  // ...done
}
