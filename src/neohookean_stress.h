#ifndef __NEOHOOKEAN_STRESS_H__
#define __NEOHOOKEAN_STRESS_H__

inline
void Neohookean::calculate_stress ()
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

  // Cauchy stress tensor
  for (unsigned int i=0; i<3; ++i) {
    for (unsigned int j=0; j<3; ++j) {
      if (i == j) {
      	this->sigma(i,i) = (mu/J) * b(i,i) - (mu/J) + (0.5*lambda*(J-1.0/J));
      } else {
      	this->sigma(i,j) = (mu/J) * b(i,j);
      }
    }
  }
}
#endif // __NEOHOOKEAN_STRESS_H__
