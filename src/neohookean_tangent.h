#ifndef __NEOHOOKEAN_TANGENT_H__
#define __NEOHOOKEAN_TANGENT_H__

inline
void Neohookean::calculate_tangent ()
{
  // Lame material parameters
  const Real mu = 0.5*Young/(1.0+Poisson);
  const Real lambda = Young*Poisson/((1.0+Poisson)*(1.0-2.0*Poisson));

  // identity tensor
  RealTensor I;
  I(0,0) = I(1,1) = I(2,2) = 1.0;

  // plastic (inelastic) deformation gradient tensor
  RealTensor Fp;
  Fp = I;
  // elastic deformation gradient tensor
  RealTensor Fe;
  Fe = F * inverse(Fp);

  // left Green-Cauchy deformation tensor
  const RealTensor b = (this->F*this->F.transpose());

  const Real J = this->F.det();
  const Real I1 = trace(b);

  this->tangent.resize(6, 6);

  // tangent stiffness tensor (in Voigt form) with respect
  // to the current (deformed) configuration
  for (unsigned int i=0; i<3; ++i) {
    for (unsigned int j=0; j<3; ++j) {
      if (i == j) {
        this->tangent(i,i) = 2.0*mu/J + lambda/J;
        this->tangent(i+3,i+3) = mu/J - 0.5*lambda*(J-1.0/J);
      } else {
        this->tangent(i,j) = lambda*J;
      }
    }
  }
}
#endif // __NEOHOOKEAN_TANGENT_H__
