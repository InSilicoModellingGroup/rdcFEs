// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Source has been modified based on the original code produced by Robert Weidlich (2012)


#include "./solid_system.h"

using namespace libMesh;

void SolidSystem::save_initial_mesh ()
{
  EquationSystems & es = this->get_equation_systems();

  System & aux_sys = es.get_system("SolidSystem::auxiliary");

  // loop over all nodes and copy the location from this system to
  // the auxiliary system
  for (const auto & node : this->get_mesh().local_node_ptr_range())
    for (unsigned int d=0; d<3; ++d)
      {
        const unsigned int source_dof = node->dof_number(this->number(), this->var[d], 0);
        const unsigned int dest_dof = node->dof_number(aux_sys.number(), this->undefo_var[d], 0);

        const Number value = this->current_local_solution->el(source_dof);
        aux_sys.current_local_solution->set(dest_dof, value);
      }
  aux_sys.current_local_solution->close();
  // fill global solution vector from local ones
  (*aux_sys.solution) = (*aux_sys.current_local_solution);
  aux_sys.solution->close();
}

void SolidSystem::init_data ()
{
  EquationSystems & es = this->get_equation_systems();

  // add the node positions as primary variables
  this->var[0] = this->variable_number("x");
  this->var[1] = this->variable_number("y");
  this->var[2] = this->variable_number("z");

  // add variables for storing the initial mesh to an auxiliary system
  ExplicitSystem & aux_sys = es.get_system<ExplicitSystem>("SolidSystem::auxiliary");
  this->undefo_var[0] = aux_sys.variable_number("undefo_x");
  this->undefo_var[1] = aux_sys.variable_number("undefo_y");
  this->undefo_var[2] = aux_sys.variable_number("undefo_z");

  // set the time stepping for the solid solver
  this->deltat = es.parameters.get<Real>("time_step");

  // so the parent's initialization after variables are defined
  FEMSystem::init_data();

  // tell the system which variables are containing the node positions
  set_mesh_system(this);

  this->set_mesh_x_var(this->var[0]);
  this->set_mesh_y_var(this->var[1]);
  this->set_mesh_z_var(this->var[2]);

  // fill the variables with the position of the nodes
  this->mesh_position_get();

  System::reinit();

  // set some options for the 'DiffSolver'
  DiffSolver & solver = *(this->time_solver->diff_solver().get());
  solver.quiet = es.parameters.get<bool>("solver/quiet");
  solver.max_nonlinear_iterations = es.parameters.get<int>("solver/nonlinear/max_nonlinear_iterations");
  solver.relative_step_tolerance = es.parameters.get<Real>("solver/nonlinear/relative_step_tolerance");
  solver.relative_residual_tolerance = es.parameters.get<Real>("solver/nonlinear/relative_residual_tolerance");
  solver.absolute_residual_tolerance = es.parameters.get<Real>("solver/nonlinear/absolute_residual_tolerance");
  solver.verbose = (! es.parameters.get<bool>("solver/quiet"));

  ((NewtonSolver &) solver).require_residual_reduction = es.parameters.get<bool>("solver/nonlinear/require_reduction");

  // set the linear solver options
  solver.max_linear_iterations = es.parameters.get<int>("solver/linear/max_linear_iterations");
  solver.initial_linear_tolerance = es.parameters.get<Real>("solver/linear/initial_linear_tolerance");
}

void SolidSystem::update ()
{
  EquationSystems & es = this->get_equation_systems();

  System::update();
  this->mesh_position_set();

  System & disp_sys = es.get_system("SolidSystem::displacement");
  System & aux_sys = es.get_system("SolidSystem::auxiliary");

  disp_sys.current_local_solution->close();
  aux_sys.current_local_solution->close();
  this->current_local_solution->close();
  //
  (*disp_sys.current_local_solution) = (*this->current_local_solution);
  disp_sys.current_local_solution->add(-1.0, *aux_sys.current_local_solution);
  disp_sys.current_local_solution->close();
  //
  (*disp_sys.solution) = (*disp_sys.current_local_solution);
  disp_sys.solution->close();
}

void SolidSystem::init_context (DiffContext& context)
{
  FEMContext& c = cast_ref<FEMContext&>(context);

  // Pre-request all the data needed
  FEBase* elem_fe = nullptr;
  c.get_element_fe(0, elem_fe);

  elem_fe->get_JxW();
  elem_fe->get_phi();
  elem_fe->get_dphi();
  elem_fe->get_xyz();

  FEBase* side_fe = nullptr;
  c.get_side_fe(0, side_fe);

  side_fe->get_JxW();
  side_fe->get_phi();
  side_fe->get_xyz();
}

bool SolidSystem::element_time_derivative (bool request_jacobian,
                                           DiffContext& context)
{
  // compute contribution from internal forces in elem_residual and contribution
  // from linearized internal forces (stiffness matrix) in 'elem_jacobian'
  EquationSystems & es = this->get_equation_systems();

  FEMContext& c = cast_ref<FEMContext&>(context);

  const Elem& elem = c.get_elem();

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase* elem_fe = nullptr;
  c.get_element_fe(0, elem_fe);

  const std::vector<Real>& JxW = elem_fe->get_JxW();
  const std::vector<std::vector<RealGradient>>& dphi = elem_fe->get_dphi();

  const unsigned int n_u_dofs = c.n_dof_indices(this->var[0]);
  libmesh_assert(n_u_dofs == c.n_dof_indices(this->var[1]));
  libmesh_assert(n_u_dofs == c.n_dof_indices(this->var[2]));

  const unsigned int n_qpoints = c.get_element_qrule().n_points();

  DenseMatrix<Real> stiff;
  DenseVector<Real> res;

  const std::string material_ID = std::to_string(elem.subdomain_id());
  const Real E = es.parameters.get<Real>("material/"+material_ID+"/Neohookean/Young"),
             v = es.parameters.get<Real>("material/"+material_ID+"/Neohookean/Poisson");
  Neohookean material(dphi, E, v);

  TransientExplicitSystem & aux_system =
    es.get_system<TransientExplicitSystem>("SolidSystem::auxiliary");

  // assume symmetry of local stiffness matrices
  const bool use_symmetry = es.parameters.get<bool>("solver/assembly_use_symmetry");

  // build the element Jacobian and residual, calculated at each
  // quadrature point by summing the solution degree-of-freedom values by
  // the appropriate weight functions.
  // This class just takes care of the assembly. The matrix of
  // the jacobian and the residual vector are provided by the
  // constitutive formulation.

  std::vector<dof_id_type> undefo_dofs[3];
  for (unsigned int d=0; d<3; ++d)
    aux_system.get_dof_map().dof_indices(&c.get_elem(), undefo_dofs[d], this->undefo_var[d]);

  for (unsigned int qp=0; qp!=n_qpoints; qp++)
    {
      // Inverse deformation gradient tensor
      VectorValue<Gradient> grad_X;
      for (unsigned int d=0; d<3; ++d)
        {
          std::vector<Number> XYZ_undefo;
          aux_system.current_local_solution->get(undefo_dofs[d], XYZ_undefo);

          for (unsigned int l=0; l!=n_u_dofs; l++)
            grad_X(d).add_scaled(dphi[l][qp], XYZ_undefo[l]);
        }

      // Initialize the constitutive formulation with the current displacement
      // gradient
      material.init_for_qp(qp, grad_X, request_jacobian);

      // Acquire, scale and assemble residual and stiffness
      for (unsigned int i=0; i<n_u_dofs; i++)
        {
          material.get_residual(res, i);
          res.scale(JxW[qp]);
          for (unsigned int ii=0; ii<3; ++ii)
            (c.get_elem_residual(ii))(i) += res(ii);

          if (request_jacobian && c.elem_solution_derivative)
            {
              libmesh_assert(1.0==c.elem_solution_derivative);

              for (unsigned int j=(use_symmetry ? i : 0); j < n_u_dofs; j++)
                {
                  material.get_linearized_stiffness(stiff, i, j);
                  stiff.scale(JxW[qp]);
                  for (unsigned int ii=0; ii<3; ++ii)
                    {
                      for (unsigned int jj=0; jj<3; ++jj)
                        {
                          (c.get_elem_jacobian(ii,jj))(i, j) += stiff(ii, jj);
                          if (use_symmetry && i != j)
                            (c.get_elem_jacobian(ii,jj))(j, i) += stiff(jj, ii);
                        }
                    }
                }
            }
        }
    } // end of the quadrature point qp-loop

  return request_jacobian;
}

bool SolidSystem::side_time_derivative (bool request_jacobian,
                                        DiffContext& context)
{
  EquationSystems & es = this->get_equation_systems();

  FEMContext& c = cast_ref<FEMContext&>(context);

  const Elem& elem = c.get_elem();

  TransientExplicitSystem & aux_system =
    es.get_system<TransientExplicitSystem>("SolidSystem::auxiliary");

  // BCs are stored in the simulation parameters as array containing sequences of
  // 4 numbers: Id of the side for the displacements and three values describing the
  // displacement, e.g.: bc/displacement = '5 nan nan -1.0'. This will move all nodes of
  // side 5 about 1.0 units down the z-axis while leaving all other directions unrestricted

  // get the current load step
  const Real ratio = es.parameters.get<Real>("time")
                   * 1.000001;

  const std::set<int> BCs_set = export_integers(es.parameters.get<std::string>("BCs"));
  // apply Dirichlet (displacement) boundary conditions using the penalty method
  for (auto bc : BCs_set)
    {
      // get ID of the side for this BC
      const boundary_id_type ID = cast_int<boundary_id_type>(bc);
      // check if current side of FE is on the boundary that needs being restricted
      if ( ! es.get_mesh().get_boundary_info().has_boundary_id(&c.get_elem(), c.get_side(), ID) )
        continue;

      Point diff_value(es.parameters.get<Point>("BC/"+std::to_string(bc)+"/displacement"));
      // ...now scale according to current load step
      diff_value *= ratio;

      const Real penalty_number = es.parameters.get<Real>("BCs/displacement_penalty");

      FEBase* fe = nullptr;
      c.get_side_fe(0, fe);

      const std::vector<std::vector<Real>>& phi = fe->get_phi();
      const std::vector<Real>& JxW = fe->get_JxW();
      const std::vector<Point>& coords = fe->get_xyz();

      const unsigned int n_x_dofs = c.n_dof_indices(this->var[0]);

      std::vector<dof_id_type> undefo_dofs[3];
      for (unsigned int d=0; d<3; ++d)
        aux_system.get_dof_map().dof_indices(&c.get_elem(), undefo_dofs[d], this->undefo_var[d]);

      for (unsigned int qp=0; qp<c.get_side_qrule().n_points(); ++qp)
        {
          // Calculate coordinates of the quadrature point on
          // undeformed mesh (reference configuration)
          Point orig_point;
          for (unsigned int i=0; i<n_x_dofs; ++i)
            {
              for (unsigned int d=0; d<3; ++d)
                {
                  Number orig_val = aux_system.current_solution(undefo_dofs[d][i]);
                  orig_point(d) += phi[i][qp] * orig_val;
                }
            }

          // Calculate displacement to be enforced.
          const Point diff = coords[qp] - orig_point
                           - diff_value;

          // Assemble
          for (unsigned int i=0; i<n_x_dofs; ++i)
            {
              for (unsigned int di=0; di<3; ++di)
                {
                  if (libmesh_isnan(diff(di))) continue;

                  const Real val = JxW[qp] * phi[i][qp] * diff(di)
                                 * penalty_number;
                  (c.get_elem_residual(this->var[di]))(i) += val;
                }
              if (request_jacobian)
                {
                  for (unsigned int j=0; j<n_x_dofs; ++j)
                    {
                      for (unsigned int dj=0; dj<3; ++dj)
                        {
                          if (libmesh_isnan(diff(dj))) continue;

                          const Real val = JxW[qp] * phi[i][qp] * phi[j][qp]
                                         * penalty_number;
                          (c.get_elem_jacobian(this->var[dj],this->var[dj]))(i, j) += val;
                        }
                    }
                }
            }
        }
    }

  return request_jacobian;
}

void SolidSystem::post_process ()
{
  EquationSystems & es = this->get_equation_systems();

  FEType fe_type = this->get_dof_map().variable_type(0);

  QGauss qrule(3, THIRD);

  std::unique_ptr<FEBase> fe(FEBase::build(3, fe_type));
  fe->attach_quadrature_rule(&qrule);

  ExplicitSystem& aux_system =
    es.get_system<TransientExplicitSystem>("SolidSystem::auxiliary");
  libmesh_assert_equal_to(aux_system.n_vars(), 3);

  ExplicitSystem& press_sys =
    es.get_system<ExplicitSystem>("SolidSystem::pressure");
  libmesh_assert_equal_to(press_sys.n_vars(), 1);

  // loop over the local active FEs in the mesh.
  for (const auto & elem : this->get_mesh().active_local_element_ptr_range())
    {
      std::vector<dof_id_type> dof_indices;
      this->get_dof_map().dof_indices(elem, dof_indices);

      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());

      std::vector<dof_id_type> dof_indices_x;
      this->get_dof_map().dof_indices(elem, dof_indices_x, this->var[0]);
      std::vector<dof_id_type> dof_indices_y;
      this->get_dof_map().dof_indices(elem, dof_indices_y, this->var[1]);
      std::vector<dof_id_type> dof_indices_z;
      this->get_dof_map().dof_indices(elem, dof_indices_z, this->var[2]);

      std::vector<dof_id_type> dof_indices_P;
      press_sys.get_dof_map().dof_indices(elem, dof_indices_P);

      const unsigned int n_u_dofs = dof_indices_x.size();
      libmesh_assert(n_u_dofs == dof_indices_y.size());
      libmesh_assert(n_u_dofs == dof_indices_z.size());

      const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();
      fe->reinit(elem);

      const std::string material_ID = std::to_string(elem->subdomain_id());
      const Real E = es.parameters.get<Real>("material/"+material_ID+"/Neohookean/Young"),
                 v = es.parameters.get<Real>("material/"+material_ID+"/Neohookean/Poisson");
      Neohookean material(dphi, E, v);

      // average Cauchy stress tensor calculated for the FE
      RealTensor stress_Cauchy;

      std::vector<dof_id_type> undefo_dofs[3];
      for (unsigned int d=0; d<3; ++d)
        aux_system.get_dof_map().dof_indices(elem, undefo_dofs[d], this->undefo_var[d]);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // inverse deformation gradient tensor
          VectorValue<Gradient> grad_X;
          for (unsigned int d=0; d<3; ++d)
            {
              std::vector<Number> XYZ_undefo;
              aux_system.current_local_solution->get(undefo_dofs[d], XYZ_undefo);

              for (unsigned int l=0; l!=n_u_dofs; l++)
                grad_X(d).add_scaled(dphi[l][qp], XYZ_undefo[l]);
            }

          // initialize the constitutive formulation with the current displacement
          // gradient
          material.init_for_qp(qp, grad_X, false);

          stress_Cauchy += material.get_stress_tensor();
        }

      // take the average Cauchy stress
      stress_Cauchy /= qrule.n_points();
      // calculate the mean solid stress
      const Real mss = (stress_Cauchy(0,0)+stress_Cauchy(1,1)+stress_Cauchy(2,2))/3.0;

      press_sys.solution->set(dof_indices_P[0], mss);
      // ...end this loop for active local FEs
    }

    press_sys.solution->close();
}

void SolidSystem::update_auxiliary ()
{
  EquationSystems & es = this->get_equation_systems();

  TransientExplicitSystem & aux_system =
    es.get_system<TransientExplicitSystem>("SolidSystem::auxiliary");

  // close all solution vector containers
  aux_system.current_local_solution->close();
  aux_system.old_local_solution->close();
  // copy the 'old' to the 'older' solution vector
  (*aux_system.older_local_solution) = (*aux_system.old_local_solution);
  // copy the 'current' to the 'old' solution vector
  (*aux_system.old_local_solution) = (*aux_system.current_local_solution);
}
