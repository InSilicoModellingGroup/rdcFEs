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


#ifndef __NEOHOOKEAN_H__
#define __NEOHOOKEAN_H__

#include "./utils.h"

//-------------------------------------------------------------------------------------------------
class Neohookean
{
public:
  // Constructor
  Neohookean (const std::vector<std::vector<RealGradient>>& dNdx, Real E, Real v, Real K =0.0) :
    dphi(dNdx), Young(E), Poisson(v), FibreStiffness(K) {}

  inline
  void init_for_qp (unsigned int qp, VectorValue<Gradient>& gradX, const Real* lambda,
                    const RealVectorValue& f, bool calc_tangent =false)
  {
    this->current_qp = qp;
    // initialize the class for the given displacement gradient
    // at the specified quadrature point
    RealTensor dX_dy;
    for (unsigned int i=0; i<3; ++i)
      for (unsigned int j=0; j<3; ++j)
        dX_dy(i, j) = gradX(i)(j);

    this->F = dX_dy.inverse();
    libmesh_assert_greater(this->F.det(), -TOLERANCE);

    // plastic (inelastic) deformation gradient tensor
    this->Fp.zero();
    for (unsigned int l=0; l<3; l++)
      this->Fp(l,l) = lambda[l];
    // elastic deformation gradient tensor
    this->Fe = this->F*this->Fp.inverse();

    this->A = FibreStiffness>0.0 ? f.unit() : RealVectorValue(0.0, 0.0, 0.0);

    if (calc_tangent)
      this->calculate_tangent();
    this->calculate_stress();
  }

  inline
  void get_residual (DenseVector<Real>& R, unsigned int i)
  {
    // return the residual vector for the current state
    R.resize(3);

    // Cauchy stress tensor in Voigt form
    DenseVector<Real> SV(6);
    SV = { this->sigma(0,0), this->sigma(1,1), this->sigma(2,2),
           this->sigma(0,1), this->sigma(1,2), this->sigma(0,2) };
    // shape functions derivatives matrix
    DenseMatrix<Real> B_L;
    this->build_b_0_mat(i, B_L);
    B_L.vector_mult(R, SV);
  }

  inline
  void get_linearized_stiffness (DenseMatrix<Real>& D, unsigned int i, unsigned int j)
  {
    // return the tangent stiffness matrix for the current state
    D.resize(3, 3);

    // geometric non-linearity contribution
    const Real G_IK = dphi[i][current_qp] * this->sigma * dphi[j][current_qp];
    D(0,0) += G_IK;
    D(1,1) += G_IK;
    D(2,2) += G_IK;
    // shape functions derivatives matrix
    DenseMatrix<Real> B_L;
    this->build_b_0_mat(i, B_L);
    // shape functions derivatives matrix
    DenseMatrix<Real> B_K;
    this->build_b_0_mat(j, B_K);
    B_L.right_multiply(this->tangent);
    B_L.right_multiply_transpose(B_K);
    // material non-linearity contribution
    D += B_L;
  }

  inline
  const RealTensor& get_deformation_gradient_tensor () const
  {
    // return the deformation gradient tensor - current configuration
    return this->F;
  }

  inline
  const RealTensor& get_stress_tensor () const
  {
    // return the Cauchy stress tensor - current configuration
    return this->sigma;
  }

private:

  void calculate_stress ();

  void calculate_tangent ();

  inline
  void build_b_0_mat (unsigned int i, DenseMatrix<Real>& B)
  {
    B.resize(3, 6);
    B(0,0) = this->dphi[i][this->current_qp](0);
    B(1,1) = this->dphi[i][this->current_qp](1);
    B(2,2) = this->dphi[i][this->current_qp](2);
    B(0,3) = this->dphi[i][this->current_qp](1);
    B(1,3) = this->dphi[i][this->current_qp](0);
    B(1,4) = this->dphi[i][this->current_qp](2);
    B(2,4) = this->dphi[i][this->current_qp](1);
    B(0,5) = this->dphi[i][this->current_qp](2);
    B(2,5) = this->dphi[i][this->current_qp](0);
  }

private:
  // parameters (material constants) of the Neo-Hookean hyperelastic constitutive model
  Real Young, Poisson, FibreStiffness;
  // fibre direction (unit) vector - reference configuration
  RealVectorValue A;
  // deformation gradient tensor
  RealTensor F;
  // elastic & plastic deformation gradient tensor
  RealTensor Fe, Fp;
  // Cauchy stress tensor
  RealTensor sigma;
  // tangent stiffness tensor - expressed in Voigt form
  DenseMatrix<Real> tangent;
  //
  unsigned int current_qp;
  const std::vector< std::vector<RealGradient> >& dphi;
};
//-------------------------------------------------------------------------------------------------

#include "./neohookean_stress.h"
#include "./neohookean_tangent.h"

#endif // __NEOHOOKEAN_H__
