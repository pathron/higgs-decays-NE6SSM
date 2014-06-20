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

// File generated at Fri 20 Jun 2014 21:37:45

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_TKappa_one_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = oneOver16PiSqr*(3*traceKappaAdjKappa*TKappa + 2*AbsSqr(
      Lambdax)*TKappa + AbsSqr(Sigmax)*TKappa - 0.26666666666666666*Sqr(g1)*
      TKappa - 0.65*Sqr(g1p)*TKappa - 5.333333333333333*Sqr(g3)*TKappa - 0.05*
      Sqr(g1p)*Sqr(QS)*TKappa + 0.03333333333333333*Kappa*(180*
      traceAdjKappaTKappa + 16*MassB*Sqr(g1) + 39*MassBp*Sqr(g1p) + 320*MassG*
      Sqr(g3) + 3*MassBp*Sqr(g1p)*Sqr(QS) + 120*Conj(Lambdax)*TLambdax + 60*
      Conj(Sigmax)*TSigmax) + 3*(Kappa*(Kappa).adjoint()*TKappa) + 3*(TKappa*(
      Kappa).adjoint()*Kappa));


   return beta_TKappa;
}

/**
 * Calculates the two-loop beta function of TKappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_TKappa_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = twoLoop*(2.2755555555555556*Power(g1,4)*TKappa + 4.28875
      *Power(g1p,4)*TKappa + 14.222222222222221*Power(g3,4)*TKappa + 0.005*
      Power(g1p,4)*Power(QS,4)*TKappa - 6*traceKappaAdjKappaKappaAdjKappa*
      TKappa - 6*traceYdAdjYd*AbsSqr(Lambdax)*TKappa - 2*traceYeAdjYe*AbsSqr(
      Lambdax)*TKappa - 6*traceYuAdjYu*AbsSqr(Lambdax)*TKappa - 2*AbsSqr(
      KappaPr)*AbsSqr(Sigmax)*TKappa + 0.8*traceKappaAdjKappa*Sqr(g1)*TKappa +
      1.2*AbsSqr(Lambdax)*Sqr(g1)*TKappa + 1.95*traceKappaAdjKappa*Sqr(g1p)*
      TKappa + 1.3*AbsSqr(Lambdax)*Sqr(g1p)*TKappa + 0.17333333333333334*Sqr(g1
      )*Sqr(g1p)*TKappa + 6*AbsSqr(Lambdax)*Sqr(g2)*TKappa + 16*
      traceKappaAdjKappa*Sqr(g3)*TKappa + 1.4222222222222223*Sqr(g1)*Sqr(g3)*
      TKappa + 3.466666666666667*Sqr(g1p)*Sqr(g3)*TKappa + 0.34375*Power(g1p,4)
      *Sqr(QS)*TKappa - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS)*TKappa - 0.1*
      AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TKappa - 4*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TKappa - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)*TKappa -
      0.0005555555555555556*Kappa*(16384*Power(g1,4)*MassB + 30879*Power(g1p,4)
      *MassBp + 102400*Power(g3,4)*MassG + 36*Power(g1p,4)*MassBp*Power(QS,4) +
      43200*traceKappaAdjKappaTKappaAdjKappa - 2880*traceAdjKappaTKappa*Sqr(g1
      ) + 2880*MassB*traceKappaAdjKappa*Sqr(g1) - 7020*traceAdjKappaTKappa*Sqr(
      g1p) + 7020*MassBp*traceKappaAdjKappa*Sqr(g1p) + 624*MassB*Sqr(g1)*Sqr(
      g1p) + 624*MassBp*Sqr(g1)*Sqr(g1p) - 57600*traceAdjKappaTKappa*Sqr(g3) +
      57600*MassG*traceKappaAdjKappa*Sqr(g3) + 5120*MassB*Sqr(g1)*Sqr(g3) +
      5120*MassG*Sqr(g1)*Sqr(g3) + 12480*MassBp*Sqr(g1p)*Sqr(g3) + 12480*MassG*
      Sqr(g1p)*Sqr(g3) + 2475*Power(g1p,4)*MassBp*Sqr(QS) + 540*
      traceAdjKappaTKappa*Sqr(g1p)*Sqr(QS) - 540*MassBp*traceKappaAdjKappa*Sqr(
      g1p)*Sqr(QS) + 28800*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 360*Conj(
      Lambdax)*(Lambdax*(60*traceAdjYdTYd + 20*traceAdjYeTYe + 60*traceAdjYuTYu
      + 12*MassB*Sqr(g1) + 13*MassBp*Sqr(g1p) + 60*MassWB*Sqr(g2) - MassBp*Sqr
      (g1p)*Sqr(QS)) + (60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu -
      12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2) + Sqr(g1p)*Sqr(QS))*TLambdax) +
      14400*Sigmax*Sqr(Conj(Sigmax))*TSigmax + 7200*Conj(KappaPr)*Conj(Sigmax)*
      (Sigmax*TKappaPr + KappaPr*TSigmax)) - 0.2*(60*traceAdjKappaTKappa +
      MassBp*Sqr(g1p)*Sqr(QS) + 40*Conj(Lambdax)*TLambdax + 20*Conj(Sigmax)*
      TSigmax)*(Kappa*(Kappa).adjoint()*Kappa) - 9*traceKappaAdjKappa*(Kappa*(
      Kappa).adjoint()*TKappa) - 6*AbsSqr(Lambdax)*(Kappa*(Kappa).adjoint()*
      TKappa) - 3*AbsSqr(Sigmax)*(Kappa*(Kappa).adjoint()*TKappa) - 0.25*Sqr(
      g1p)*(Kappa*(Kappa).adjoint()*TKappa) + 0.15*Sqr(g1p)*Sqr(QS)*(Kappa*(
      Kappa).adjoint()*TKappa) - 9*traceKappaAdjKappa*(TKappa*(Kappa).adjoint()
      *Kappa) - 6*AbsSqr(Lambdax)*(TKappa*(Kappa).adjoint()*Kappa) - 3*AbsSqr(
      Sigmax)*(TKappa*(Kappa).adjoint()*Kappa) + 0.25*Sqr(g1p)*(TKappa*(Kappa)
      .adjoint()*Kappa) + 0.15*Sqr(g1p)*Sqr(QS)*(TKappa*(Kappa).adjoint()*Kappa
      ) - 3*(Kappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*TKappa) - 4*(Kappa
      *(Kappa).adjoint()*TKappa*(Kappa).adjoint()*Kappa) - 3*(TKappa*(Kappa)
      .adjoint()*Kappa*(Kappa).adjoint()*Kappa));


   return beta_TKappa;
}

} // namespace flexiblesusy
