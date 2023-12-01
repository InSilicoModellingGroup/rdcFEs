#ifndef __NEOHOOKEAN_STRESS_H__
#define __NEOHOOKEAN_STRESS_H__

inline
void Neohookean::calculate_stress ()
{
  // Lame material parameters
  const Real mu = 0.5*Young/(1.0+Poisson);
  const Real lambda = Young*Poisson/((1.0+Poisson)*(1.0-2.0*Poisson));
  // fibre (axial stiffness) material parameter
  const Real koppa = FibreStiffness / 2.0;

  // right Green-Cauchy deformation tensor - elastic
  const RealTensor Ce = (this->Fe.transpose()*this->Fe);
  //
  const RealTensor CeINV = Ce.inverse();
  const RealTensor delta = RealTensor(1,0,0,0,1,0,0,0,1);
  //
  const Real I1 = Ce.tr();
  const Real I2 = 0.5*(pow2(I1)-(Ce*Ce).tr());
  const Real Je = this->Fe.det();
  const Real I4 = (this->A*Ce*this->A);
  //
  const Real J_recip = 1.0/this->F.det();
  //
  const Real dWdI1 = (mu/2.0);
  const Real dWdI2 = 0.0;
  const Real dWdJe = (-mu/Je) + (lambda/2.0*Je-lambda/2.0/Je);
  const Real dWdI4 = 0.0;
  //
  Real dI1dCe[3][3];
  Real dI2dCe[3][3];
  Real dJedCe[3][3];
  Real dI4dCe[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      dI1dCe[i][j] = delta(i,j);
      dI2dCe[i][j] = delta(i,j)*I1-Ce(i,j);
      dJedCe[i][j] = 0.5*Je*CeINV(i,j);
      dI4dCe[i][j] = A(i)*A(j);
    }
  }
  //
  Real S2pk[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      S2pk[i][j] = 2.0*dWdI1*dI1dCe[i][j]
                  +2.0*dWdI2*dI2dCe[i][j]
                  +2.0*dWdJe*dJedCe[i][j]
                  +2.0*dWdI4*dI4dCe[i][j];
    }
  }
  //
  Real Sc[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      Real F_F_S2pk(0.0);
      for (int I=0; I<3; I++)
        for (int J=0; J<3; J++)
          F_F_S2pk += F(i,I)*F(j,J)*S2pk[I][J];
      Sc[i][j] = F_F_S2pk*J_recip;
    }
  }
  //
  for (int I=0; I<3; I++)
    for (int J=0; J<3; J++)
      this->sigma(I,J) = Sc[I][J];
}
#endif // __NEOHOOKEAN_STRESS_H__
