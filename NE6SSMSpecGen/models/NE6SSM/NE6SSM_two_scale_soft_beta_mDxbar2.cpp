// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Fri 13 Jun 2014 15:36:46

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mDxbar2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_mDxbar2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = oneOver16PiSqr*(2*ms2*((Kappa).transpose()*
      Kappa.conjugate()) + 2*((TKappa).transpose()*TKappa.conjugate()) +
      mDxbar2*(Kappa).transpose()*Kappa.conjugate() + 2*((Kappa).transpose()*
      mDx2*Kappa.conjugate()) + (Kappa).transpose()*Kappa.conjugate()*mDxbar2 +
      0.5163977794943222*g1*Tr11*UNITMATRIX(3) - 0.9486832980505138*g1p*Tr14*
      UNITMATRIX(3) - 0.5333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) -
      1.8*AbsSqr(MassBp)*Sqr(g1p)*UNITMATRIX(3) - 10.666666666666666*AbsSqr(
      MassG)*Sqr(g3)*UNITMATRIX(3));


   return beta_mDxbar2;
}

/**
 * Calculates the two-loop beta function of mDxbar2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_mDxbar2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = twoLoop*(-6*traceconjTKappaTpTKappa*((Kappa).transpose(
      )*Kappa.conjugate()) - 12*ms2*traceKappaAdjKappa*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceKappaAdjKappaconjmDx2*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceKappaconjmDxbar2AdjKappa*((Kappa).transpose()
      *Kappa.conjugate()) - 4*mHd2*AbsSqr(Lambdax)*((Kappa).transpose()*
      Kappa.conjugate()) - 4*mHu2*AbsSqr(Lambdax)*((Kappa).transpose()*
      Kappa.conjugate()) - 8*ms2*AbsSqr(Lambdax)*((Kappa).transpose()*
      Kappa.conjugate()) - 2*mphi2*AbsSqr(Sigmax)*((Kappa).transpose()*
      Kappa.conjugate()) - 4*ms2*AbsSqr(Sigmax)*((Kappa).transpose()*
      Kappa.conjugate()) - 2*msbar2*AbsSqr(Sigmax)*((Kappa).transpose()*
      Kappa.conjugate()) - 4*AbsSqr(TLambdax)*((Kappa).transpose()*
      Kappa.conjugate()) - 2*AbsSqr(TSigmax)*((Kappa).transpose()*
      Kappa.conjugate()) - 0.5*ms2*Sqr(g1p)*((Kappa).transpose()*
      Kappa.conjugate()) + 0.1*ms2*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*
      Kappa.conjugate()) - 6*traceAdjKappaTKappa*((Kappa).transpose()*
      TKappa.conjugate()) + 0.5*MassBp*Sqr(g1p)*((Kappa).transpose()*
      TKappa.conjugate()) - 0.1*MassBp*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*
      TKappa.conjugate()) - 4*Conj(Lambdax)*TLambdax*((Kappa).transpose()*
      TKappa.conjugate()) - 2*Conj(Sigmax)*TSigmax*((Kappa).transpose()*
      TKappa.conjugate()) - 6*traceconjTKappaTpKappa*((TKappa).transpose()*
      Kappa.conjugate()) - 4*Conj(TLambdax)*Lambdax*((TKappa).transpose()*
      Kappa.conjugate()) - 2*Conj(TSigmax)*Sigmax*((TKappa).transpose()*
      Kappa.conjugate()) - 6*traceKappaAdjKappa*((TKappa).transpose()*
      TKappa.conjugate()) - 4*AbsSqr(Lambdax)*((TKappa).transpose()*
      TKappa.conjugate()) - 2*AbsSqr(Sigmax)*((TKappa).transpose()*
      TKappa.conjugate()) - 0.5*Sqr(g1p)*((TKappa).transpose()*TKappa.conjugate
      ()) + 0.1*Sqr(g1p)*Sqr(QS)*((TKappa).transpose()*TKappa.conjugate()) - 3*
      traceKappaAdjKappa*(mDxbar2*(Kappa).transpose()*Kappa.conjugate()) - 2*
      AbsSqr(Lambdax)*(mDxbar2*(Kappa).transpose()*Kappa.conjugate()) - AbsSqr(
      Sigmax)*(mDxbar2*(Kappa).transpose()*Kappa.conjugate()) - 0.25*Sqr(g1p)*(
      mDxbar2*(Kappa).transpose()*Kappa.conjugate()) + 0.05*Sqr(g1p)*Sqr(QS)*(
      mDxbar2*(Kappa).transpose()*Kappa.conjugate()) - 6*traceKappaAdjKappa*((
      Kappa).transpose()*mDx2*Kappa.conjugate()) - 4*AbsSqr(Lambdax)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 2*AbsSqr(Sigmax)*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - 0.5*Sqr(g1p)*((Kappa).transpose()*
      mDx2*Kappa.conjugate()) + 0.1*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*mDx2*
      Kappa.conjugate()) - 3*traceKappaAdjKappa*((Kappa).transpose()*
      Kappa.conjugate()*mDxbar2) - 2*AbsSqr(Lambdax)*((Kappa).transpose()*
      Kappa.conjugate()*mDxbar2) - AbsSqr(Sigmax)*((Kappa).transpose()*
      Kappa.conjugate()*mDxbar2) - 0.25*Sqr(g1p)*((Kappa).transpose()*
      Kappa.conjugate()*mDxbar2) + 0.05*Sqr(g1p)*Sqr(QS)*((Kappa).transpose()*
      Kappa.conjugate()*mDxbar2) - 4*ms2*((Kappa).transpose()*Kappa.conjugate()
      *(Kappa).transpose()*Kappa.conjugate()) - 2*((Kappa).transpose()*
      Kappa.conjugate()*(TKappa).transpose()*TKappa.conjugate()) - 2*((Kappa)
      .transpose()*TKappa.conjugate()*(TKappa).transpose()*Kappa.conjugate()) -
      2*((TKappa).transpose()*Kappa.conjugate()*(Kappa).transpose()*
      TKappa.conjugate()) - 2*((TKappa).transpose()*TKappa.conjugate()*(Kappa)
      .transpose()*Kappa.conjugate()) - mDxbar2*(Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate() - 2*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()
      ) - 2*((Kappa).transpose()*Kappa.conjugate()*mDxbar2*(Kappa).transpose()*
      Kappa.conjugate()) - 2*((Kappa).transpose()*Kappa.conjugate()*(Kappa)
      .transpose()*mDx2*Kappa.conjugate()) - (Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()*mDxbar2 +
      10.666666666666666*Power(g3,4)*Tr23*UNITMATRIX(3) - 0.9797958971132712*g1
      *g1p*Tr2U114*UNITMATRIX(3) - 0.9797958971132712*g1*g1p*Tr2U141*UNITMATRIX
      (3) + 2.065591117977289*g1*Tr31*UNITMATRIX(3) - 3.794733192202055*g1p*
      Tr34*UNITMATRIX(3) + 53.333333333333336*Power(g3,4)*AbsSqr(MassG)*
      UNITMATRIX(3) + 0.5333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 1.8*
      Tr2U144*Sqr(g1p)*UNITMATRIX(3) + 2.8444444444444446*AbsSqr(MassG)*Sqr(g1)
      *Sqr(g3)*UNITMATRIX(3) + 1.4222222222222223*MassB*Conj(MassG)*Sqr(g1)*Sqr
      (g3)*UNITMATRIX(3) + 9.6*AbsSqr(MassG)*Sqr(g1p)*Sqr(g3)*UNITMATRIX(3) +
      4.8*MassBp*Conj(MassG)*Sqr(g1p)*Sqr(g3)*UNITMATRIX(3) +
      0.035555555555555556*Conj(MassB)*Sqr(g1)*(384*MassB*Sqr(g1) + 27*(2*MassB
      + MassBp)*Sqr(g1p) + 40*(2*MassB + MassG)*Sqr(g3))*UNITMATRIX(3) + 0.005
      *Conj(MassBp)*Sqr(g1p)*(20*(-5 + Sqr(QS))*(2*MassBp*((Kappa).transpose()*
      Kappa.conjugate()) - (TKappa).transpose()*Kappa.conjugate()) + 3*(64*(
      MassB + 2*MassBp)*Sqr(g1) + 320*(2*MassBp + MassG)*Sqr(g3) + 9*MassBp*Sqr
      (g1p)*(267 + 2*Sqr(QS)))*UNITMATRIX(3)));


   return beta_mDxbar2;
}

} // namespace flexiblesusy
