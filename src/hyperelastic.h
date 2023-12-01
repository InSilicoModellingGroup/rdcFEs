#ifndef __NEOHOOKEAN_H__
#define __NEOHOOKEAN_H__

#include "./utils.h"

//-------------------------------------------------------------------------------------------------
class Hyperelastic
{
public:
  // Constructor
  Hyperelastic (const std::vector<std::vector<RealGradient>>& dNdx, Real E, Real v, Real K =0.0) :
    dphi(dNdx), Young(E), Poisson(v), FibreStiffness(K)
  {
    this->Voigt_map.resize(6, 2);
    //
    this->Voigt_map(0,0) = 0; this->Voigt_map(0,1) = 0;
    this->Voigt_map(1,0) = 1; this->Voigt_map(1,1) = 1;
    this->Voigt_map(2,0) = 2; this->Voigt_map(2,1) = 2;
    this->Voigt_map(3,0) = 0; this->Voigt_map(3,1) = 1;
    this->Voigt_map(4,0) = 1; this->Voigt_map(4,1) = 2;
    this->Voigt_map(5,0) = 0; this->Voigt_map(5,1) = 2;
  }

  inline
  void initialize (unsigned int qp, VectorValue<Gradient>& gradX, const Real* lambda,
                   const RealVectorValue& f, bool calculate_tangent =false)
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

    this->calculate_stress(calculate_tangent);
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
    this->build_B_matrix(i, B_L);
    B_L.vector_mult(R, SV);
  }

  inline
  void get_linearized_stiffness (DenseMatrix<Real>& D, unsigned int i, unsigned int j)
  {
    // return the tangent stiffness matrix for the current state
    D.resize(3, 3);

    // geometric non-linearity contribution
    const Real G_NN = dphi[i][current_qp] * this->sigma * dphi[j][current_qp];
    for (unsigned int n=0; n<3; n++)
      D(n,n) += G_NN;
    // shape functions derivatives matrix
    DenseMatrix<Real> B_L;
    this->build_B_matrix(i, B_L);
    // shape functions derivatives matrix
    DenseMatrix<Real> B_K;
    this->build_B_matrix(j, B_K);
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

  void build_B_matrix (unsigned int i, DenseMatrix<Real>& B);

  void calculate_stress (bool calculate_tangent);

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
  // map of indices into Voigt form
  DenseMatrix<Real> Voigt_map;
  //
  unsigned int current_qp;
  const std::vector< std::vector<RealGradient> >& dphi;
};
//-------------------------------------------------------------------------------------------------

#include "./hyperlastic_inline.h"

#endif // __NEOHOOKEAN_H__
