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
  System& aux_sys = this->get_equation_systems().get_system("auxiliary");

  // Loop over all nodes and copy the location from the current system to
  // the auxiliary system.
  for (const auto & node : this->get_mesh().local_node_ptr_range())
    for (unsigned int d=0; d<3; ++d)
      {
        const unsigned int source_dof = node->dof_number(this->number(), this->var[d], 0);
        const unsigned int dest_dof = node->dof_number(aux_sys.number(), this->undefo_var[d], 0);

        const Number value = this->current_local_solution->el(source_dof);
        aux_sys.current_local_solution->set(dest_dof, value);
      }
}

void SolidSystem::init_data ()
{
  // get the default order of the used elements;
  // assumption: just one type of elements in the mesh
  const Order order = (*(this->get_mesh().elements_begin()))->default_order();

  // Add the node positions as primary variables.
  this->var[0] = this->add_variable("x", order);
  this->var[1] = this->add_variable("y", order);
  this->var[2] = this->add_variable("z", order);

  // Add variables for storing the initial mesh to an auxiliary system.
  System & aux_sys = this->get_equation_systems().get_system("auxiliary");
  this->undefo_var[0] = aux_sys.add_variable("undefo_x", order);
  this->undefo_var[1] = aux_sys.add_variable("undefo_y", order);
  this->undefo_var[2] = aux_sys.add_variable("undefo_z", order);

  // Add variables for storing the initial mesh to an auxiliary system.
  System & disp_sys = this->get_equation_systems().get_system("displacement");
  disp_sys.add_variable("U_x", order);
  disp_sys.add_variable("U_y", order);
  disp_sys.add_variable("U_z", order);

  // add variables for storing the initial mesh to an auxiliary system.
  System & press_sys = this->get_equation_systems().get_system("pressure");
  press_sys.add_variable("p", CONSTANT, MONOMIAL);

  // set the time stepping options
  this->deltat = this->args("solver/time_step", 0.1);

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
  solver.quiet = this->args("solver/quiet", false);
  solver.max_nonlinear_iterations = this->args("solver/nonlinear/max_nonlinear_iterations", 100);
  solver.relative_step_tolerance = this->args("solver/nonlinear/relative_step_tolerance", 1.e-3);
  solver.relative_residual_tolerance = this->args("solver/nonlinear/relative_residual_tolerance", 1.e-8);
  solver.absolute_residual_tolerance = this->args("solver/nonlinear/absolute_residual_tolerance", 1.e-8);
  solver.verbose = ! this->args("solver/quiet", false);

  ((NewtonSolver&) solver).require_residual_reduction = this->args("solver/nonlinear/require_reduction", false);

  // And the linear solver options
  solver.max_linear_iterations = this->args("max_linear_iterations", 50000);
  solver.initial_linear_tolerance = this->args("initial_linear_tolerance", 1.e-3);
}

void SolidSystem::update ()
{
  System::update();
  this->mesh_position_set();

  System & disp_sys = this->get_equation_systems().get_system("displacement");
  System & aux_sys = this->get_equation_systems().get_system("auxiliary");

  disp_sys.current_local_solution->close();
  aux_sys.current_local_solution->close();
  this->current_local_solution->close();

  (*disp_sys.current_local_solution) = (*this->current_local_solution);
  disp_sys.current_local_solution->add(-1.0, *aux_sys.current_local_solution);
  disp_sys.current_local_solution->close();

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
  FEMContext& c = cast_ref<FEMContext&>(context);

  // Access directly the finite element
  const Elem& elem = c.get_elem();

  // First we get some references to cell-specific data that
  // will be used to assemble the linear system.
  FEBase* elem_fe = nullptr;
  c.get_element_fe(0, elem_fe);

  // Element Jacobian * quadrature weights for interior integration
  const std::vector<Real>& JxW = elem_fe->get_JxW();

  // Element basis functions
  const std::vector<std::vector<RealGradient>>& dphi = elem_fe->get_dphi();

  // The number of local degrees of freedom in each variable
  const unsigned int n_u_dofs = c.n_dof_indices(this->var[0]);
  libmesh_assert(n_u_dofs == c.n_dof_indices(this->var[1]));
  libmesh_assert(n_u_dofs == c.n_dof_indices(this->var[2]));

  const unsigned int n_qpoints = c.get_element_qrule().n_points();

  // Some matrices and vectors for storing the results of the constitutive
  // law
  DenseMatrix<Real> stiff;
  DenseVector<Real> res;

  // Instantiate the constitutive law
  // Just calculate jacobian contribution when we need to
  Neohookean material(dphi, elem.subdomain_id(), this->args);

  // Get a reference to the auxiliary system
  TransientExplicitSystem& aux_system =
    this->get_equation_systems().get_system<TransientExplicitSystem>("auxiliary");

  // Assume symmetry of local stiffness matrices
  const bool use_symmetry = this->args("assembly/use_symmetry", false);

  // Now we will build the element Jacobian and residual, calculated at each
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
  FEMContext& c = cast_ref<FEMContext &>(context);

  const Elem& elem = c.get_elem();

  TransientExplicitSystem& aux_system =
    this->get_equation_systems().get_system<TransientExplicitSystem>("auxiliary");

  // The BC are stored in the simulation parameters as array containing sequences of
  // four numbers: Id of the side for the displacements and three values describing the
  // displacement. E.g.: bc/displacement = '5 nan nan -1.0'. This will move all nodes of
  // side 5 about 1.0 units down the z-axis while leaving all other directions unrestricted

  // get the current load step
  const Real ratio = this->get_equation_systems().parameters.get<Real>("progress")
                   * 1.000001;
  // get number of BCs to enforce
  boundary_id_type num_bc =
    cast_int<boundary_id_type>(this->args.vector_variable_size("bc/displacement"));

  libmesh_error_msg_if(num_bc % 4 != 0,
                       "ERROR, Odd number of values in displacement boundary condition.");
  num_bc /= 4;
  // apply Dirichlet (displacement) boundary conditions using the penalty method

  for (boundary_id_type nbc=0; nbc<num_bc; nbc++)
    {
      // get IDs of the side for this BC
      boundary_id_type positive_boundary_id =
        cast_int<boundary_id_type>(this->args("bc/displacement", 1, nbc*4));

      // current side may not be on the boundary to be restricted
      if ( ! this->get_mesh().get_boundary_info().has_boundary_id( &c.get_elem(),
                                                                    c.get_side(),
                                                                   positive_boundary_id ) )
        continue;

      Point diff_value;
      for (unsigned int d=0; d<3; ++d)
        diff_value(d) = this->args("bc/displacement", NAN, ((nbc*4+1)+d));
      // ...now scale according to current load step
      diff_value *= ratio;

      const Real penalty_number = this->args("bc/displacement_penalty", 1.0e+7);

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
  FEType fe_type = this->get_dof_map().variable_type(0);

  QGauss qrule (3, THIRD);

  std::unique_ptr<FEBase> fe(FEBase::build(3, fe_type));
  fe->attach_quadrature_rule(&qrule);

  ExplicitSystem& aux_system =
    this->get_equation_systems().get_system<TransientExplicitSystem>("auxiliary");
  libmesh_assert_equal_to(aux_system.n_vars(), 3);

  ExplicitSystem& press_sys =
    this->get_equation_systems().get_system<ExplicitSystem>("pressure");
  libmesh_assert_equal_to(aux_system.n_vars(), 1);

  // Loop over the local active FEs in the mesh.
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

      // The number of local degrees of freedom in each variable
      const unsigned int n_u_dofs = dof_indices_x.size();
      libmesh_assert(n_u_dofs == dof_indices_y.size());
      libmesh_assert(n_u_dofs == dof_indices_z.size());

      const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();
      fe->reinit (elem);

      Neohookean material(dphi, elem->subdomain_id(), this->args);

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
  TransientExplicitSystem& aux_system =
    this->get_equation_systems().get_system<TransientExplicitSystem>("auxiliary");

  aux_system.current_local_solution->close();
  aux_system.old_local_solution->close();

  // Copy the 'old' to the 'older' solution vector
  (*aux_system.older_local_solution) = (*aux_system.old_local_solution);
  // Now copy the 'current' to the 'old' solution vector
  (*aux_system.old_local_solution) = (*aux_system.current_local_solution);
}
