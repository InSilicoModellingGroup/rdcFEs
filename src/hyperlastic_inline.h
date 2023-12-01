
inline
void Hyperelastic::build_B_matrix (unsigned int i, DenseMatrix<Real>& B)
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

inline
void Hyperelastic::calculate_stress (bool calculate_tangent)
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
  const RealTensor FpINV = this->Fp.inverse();
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
  const Real dWdI4 = (-koppa);
  //
  const Real d2WdI1dI1 = 0.0;
  const Real d2WdI2dI2 = 0.0;
  const Real d2WdJedJe = (mu/Je/Je) + (lambda/2.0+lambda/2.0/Je/Je);
  const Real d2WdI4dI4 = 0.0;
  //
  Real dI1dCe[3][3];
  Real dI2dCe[3][3];
  Real dJedCe[3][3];
  Real dI4dCe[3][3];
  Real d2I2dCe2[3][3][3][3];
  Real d2JedCe2[3][3][3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      dI1dCe[i][j] = delta(i,j);
      dI2dCe[i][j] = delta(i,j)*I1-Ce(i,j);
      dJedCe[i][j] = 0.5*Je*CeINV(i,j);
      dI4dCe[i][j] = A(i)*A(j);
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          d2I2dCe2[i][j][k][l] =
            delta(i,j)*delta(k,l)-0.5*delta(i,k)*delta(j,l)-0.5*delta(i,l)*delta(j,k);
          d2JedCe2[i][j][k][l] =
            0.25*Je*CeINV(i,j)*CeINV(k,l)-0.25*Je*CeINV(i,k)*CeINV(j,l)-0.25*Je*CeINV(i,l)*CeINV(j,k);
        }
      }
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
  //
  if (!calculate_tangent) return;
  //
  Real dSdCe[3][3][3][3];
  Real dCedC[3][3][3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          dSdCe[i][j][k][l] = 4.0*dWdI2*d2I2dCe2[i][j][k][l]
                             +4.0*dWdJe*d2JedCe2[i][j][k][l]
                             +4.0*d2WdI1dI1*dI1dCe[i][j]*dI1dCe[k][l]
                             +4.0*d2WdI2dI2*dI2dCe[i][j]*dI2dCe[k][l]
                             +4.0*d2WdJedJe*dJedCe[i][j]*dJedCe[k][l]
                             +4.0*d2WdI4dI4*dI4dCe[i][j]*dI4dCe[k][l];
          dCedC[i][j][k][l] = 0.5*FpINV(k,i)*FpINV(j,l)+0.5*FpINV(l,i)*FpINV(k,j);
        }
      }
    }
  }
  //
  Real dSdC[3][3][3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          Real dSdCe__dCedC(0.0);
          for (int m=0; m<3; m++)
            for (int n=0; n<3; n++)
              dSdCe__dCedC += dSdCe[i][j][m][n]*dCedC[m][n][k][l];
          dSdC[i][j][k][l] = dSdCe__dCedC;
        }
      }
    }
  }
  //
  Real tsm[3][3][3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        for (int l=0; l<3; l++) {
          Real F_F_F_F_dSdC(0.0);
          for (int I=0; I<3; I++)
            for (int J=0; J<3; J++)
              for (int K=0; K<3; K++)
                for (int L=0; L<3; L++)
                  F_F_F_F_dSdC +=
                    F(i,I)*F(j,J)*F(k,K)*F(l,L)*dSdC[I][J][K][L];
          tsm[i][j][k][l] = F_F_F_F_dSdC*J_recip;
        }
      }
    }
  }
  // tangent stiffness matrix - current configuration
  this->tangent.resize(6, 6);
  // ...in condensed (Voigt) form
  this->tangent(0,0) = tsm[0][0][0][0];
  this->tangent(0,1) = tsm[0][0][1][1];
  this->tangent(0,2) = tsm[0][0][2][2];
  this->tangent(0,3) = tsm[0][0][0][1];
  this->tangent(0,4) = tsm[0][0][1][2];
  this->tangent(0,5) = tsm[0][0][0][2];
  this->tangent(1,0) = tsm[1][1][0][0];
  this->tangent(1,1) = tsm[1][1][1][1];
  this->tangent(1,2) = tsm[1][1][2][2];
  this->tangent(1,3) = tsm[1][1][0][1];
  this->tangent(1,4) = tsm[1][1][1][2];
  this->tangent(1,5) = tsm[1][1][0][2];
  this->tangent(2,0) = tsm[2][2][0][0];
  this->tangent(2,1) = tsm[2][2][1][1];
  this->tangent(2,2) = tsm[2][2][2][2];
  this->tangent(2,3) = tsm[2][2][0][1];
  this->tangent(2,4) = tsm[2][2][1][2];
  this->tangent(2,5) = tsm[2][2][0][2];
  this->tangent(3,0) = tsm[0][1][0][0];
  this->tangent(3,1) = tsm[0][1][1][1];
  this->tangent(3,2) = tsm[0][1][2][2];
  this->tangent(3,3) = tsm[0][1][0][1];
  this->tangent(3,4) = tsm[0][1][1][2];
  this->tangent(3,5) = tsm[0][1][0][2];
  this->tangent(4,0) = tsm[1][2][0][0];
  this->tangent(4,1) = tsm[1][2][1][1];
  this->tangent(4,2) = tsm[1][2][2][2];
  this->tangent(4,3) = tsm[1][2][0][1];
  this->tangent(4,4) = tsm[1][2][1][2];
  this->tangent(4,5) = tsm[1][2][0][2];
  this->tangent(5,0) = tsm[0][2][0][0];
  this->tangent(5,1) = tsm[0][2][1][1];
  this->tangent(5,2) = tsm[0][2][2][2];
  this->tangent(5,3) = tsm[0][2][0][1];
  this->tangent(5,4) = tsm[0][2][1][2];
  this->tangent(5,5) = tsm[0][2][0][2];
}
