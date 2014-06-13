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

// File generated at Fri 13 Jun 2014 15:36:45

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = oneOver16PiSqr*(0.31622776601683794*g1p*QS*Tr14 + 6*
      traceconjTKappaTpTKappa + 6*ms2*traceKappaAdjKappa + 6*
      traceKappaAdjKappaconjmDx2 + 6*traceKappaconjmDxbar2AdjKappa + 4*(mHd2 +
      mHu2 + ms2)*AbsSqr(Lambdax) + 2*mphi2*AbsSqr(Sigmax) + 2*ms2*AbsSqr(
      Sigmax) + 2*msbar2*AbsSqr(Sigmax) + 4*AbsSqr(TLambdax) + 2*AbsSqr(TSigmax
      ) - 0.2*AbsSqr(MassBp)*Sqr(g1p)*Sqr(QS));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjTKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa;
   const double traceKappaAdjTKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa;
   const double traceKappaAdjKappaKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2;
   const double traceKappaAdjKappaKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaconjmDx2KappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa;
   const double traceKappaconjmDxbar2AdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   beta_ms2 = 0.005*twoLoop*(Conj(MassBp)*Sqr(g1p)*(3*(20*
      traceAdjKappaTKappa*(-13 + Sqr(QS)) - 40*MassBp*traceKappaAdjKappa*(-13 +
      Sqr(QS)) + MassBp*Sqr(g1p)*Sqr(QS)*(249 + 4*Sqr(QS))) - 40*Conj(Lambdax)
      *(-13 + Sqr(QS))*(2*MassBp*Lambdax - TLambdax)) + 20*(12.649110640673518*
      g1p*QS*Tr34 - 120*traceKappaAdjKappaconjmDx2KappaAdjKappa - 240*ms2*
      traceKappaAdjKappaKappaAdjKappa - 120*
      traceKappaAdjKappaKappaAdjKappaconjmDx2 - 120*
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa - 240*
      traceKappaAdjKappaTKappaAdjTKappa - 240*traceKappaAdjTKappaTKappaAdjKappa
      - 120*traceKappaconjmDxbar2AdjKappaKappaAdjKappa - 40*AbsSqr(Sigmax)*
      AbsSqr(TKappaPr) - 120*traceYdAdjYd*AbsSqr(TLambdax) - 40*traceYeAdjYe*
      AbsSqr(TLambdax) - 120*traceYuAdjYu*AbsSqr(TLambdax) - 160*AbsSqr(Sigmax)
      *AbsSqr(TSigmax) - 120*traceAdjYdTYd*Conj(TLambdax)*Lambdax - 40*
      traceAdjYeTYe*Conj(TLambdax)*Lambdax - 120*traceAdjYuTYu*Conj(TLambdax)*
      Lambdax - 16*MassB*traceconjTKappaTpKappa*Sqr(g1) + 16*
      traceconjTKappaTpTKappa*Sqr(g1) + 16*ms2*traceKappaAdjKappa*Sqr(g1) + 16*
      traceKappaAdjKappaconjmDx2*Sqr(g1) + 16*traceKappaconjmDxbar2AdjKappa*Sqr
      (g1) + 32*traceKappaAdjKappa*AbsSqr(MassB)*Sqr(g1) + 24*AbsSqr(TLambdax)*
      Sqr(g1) - 16*traceAdjKappaTKappa*Conj(MassB)*Sqr(g1) - 24*MassB*Conj(
      TLambdax)*Lambdax*Sqr(g1) - 39*MassBp*traceconjTKappaTpKappa*Sqr(g1p) +
      39*traceconjTKappaTpTKappa*Sqr(g1p) + 39*ms2*traceKappaAdjKappa*Sqr(g1p)
      + 39*traceKappaAdjKappaconjmDx2*Sqr(g1p) + 39*
      traceKappaconjmDxbar2AdjKappa*Sqr(g1p) + 26*AbsSqr(TLambdax)*Sqr(g1p) -
      26*MassBp*Conj(TLambdax)*Lambdax*Sqr(g1p) + 120*AbsSqr(TLambdax)*Sqr(g2)
      - 120*MassWB*Conj(TLambdax)*Lambdax*Sqr(g2) - 320*MassG*
      traceconjTKappaTpKappa*Sqr(g3) + 320*traceconjTKappaTpTKappa*Sqr(g3) +
      320*ms2*traceKappaAdjKappa*Sqr(g3) + 320*traceKappaAdjKappaconjmDx2*Sqr(
      g3) + 320*traceKappaconjmDxbar2AdjKappa*Sqr(g3) + 640*traceKappaAdjKappa*
      AbsSqr(MassG)*Sqr(g3) - 320*traceAdjKappaTKappa*Conj(MassG)*Sqr(g3) + 2*
      Tr2U144*Sqr(g1p)*Sqr(QS) + 3*MassBp*traceconjTKappaTpKappa*Sqr(g1p)*Sqr(
      QS) - 3*traceconjTKappaTpTKappa*Sqr(g1p)*Sqr(QS) - 3*ms2*
      traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) - 3*traceKappaAdjKappaconjmDx2*Sqr(
      g1p)*Sqr(QS) - 3*traceKappaconjmDxbar2AdjKappa*Sqr(g1p)*Sqr(QS) - 2*
      AbsSqr(TLambdax)*Sqr(g1p)*Sqr(QS) + 2*MassBp*Conj(TLambdax)*Lambdax*Sqr(
      g1p)*Sqr(QS) - 160*(mHd2 + mHu2 + ms2)*Sqr(Conj(Lambdax))*Sqr(Lambdax) -
      80*mphi2*Sqr(Conj(Sigmax))*Sqr(Sigmax) - 80*ms2*Sqr(Conj(Sigmax))*Sqr(
      Sigmax) - 80*msbar2*Sqr(Conj(Sigmax))*Sqr(Sigmax) + 2*Conj(Lambdax)*(-60*
      traceconjTYdTpTYd*Lambdax - 20*traceconjTYeTpTYe*Lambdax - 60*
      traceconjTYuTpTYu*Lambdax - 60*tracemd2YdAdjYd*Lambdax - 20*
      traceme2YeAdjYe*Lambdax - 20*traceml2AdjYeYe*Lambdax - 60*tracemq2AdjYdYd
      *Lambdax - 60*tracemq2AdjYuYu*Lambdax - 60*tracemu2YuAdjYu*Lambdax - 120*
      mHd2*traceYdAdjYd*Lambdax - 60*mHu2*traceYdAdjYd*Lambdax - 60*ms2*
      traceYdAdjYd*Lambdax - 40*mHd2*traceYeAdjYe*Lambdax - 20*mHu2*
      traceYeAdjYe*Lambdax - 20*ms2*traceYeAdjYe*Lambdax - 60*mHd2*traceYuAdjYu
      *Lambdax - 120*mHu2*traceYuAdjYu*Lambdax - 60*ms2*traceYuAdjYu*Lambdax -
      160*AbsSqr(TLambdax)*Lambdax + 12*mHd2*Lambdax*Sqr(g1) + 12*mHu2*Lambdax*
      Sqr(g1) + 12*ms2*Lambdax*Sqr(g1) + 13*mHd2*Lambdax*Sqr(g1p) + 13*mHu2*
      Lambdax*Sqr(g1p) + 13*ms2*Lambdax*Sqr(g1p) + 60*mHd2*Lambdax*Sqr(g2) + 60
      *mHu2*Lambdax*Sqr(g2) + 60*ms2*Lambdax*Sqr(g2) - mHd2*Lambdax*Sqr(g1p)*
      Sqr(QS) - mHu2*Lambdax*Sqr(g1p)*Sqr(QS) - ms2*Lambdax*Sqr(g1p)*Sqr(QS) +
      12*Conj(MassB)*Sqr(g1)*(2*MassB*Lambdax - TLambdax) + 60*Conj(MassWB)*Sqr
      (g2)*(2*MassWB*Lambdax - TLambdax) - 60*traceconjTYdTpYd*TLambdax - 20*
      traceconjTYeTpYe*TLambdax - 60*traceconjTYuTpYu*TLambdax) - 40*Conj(
      Sigmax)*Conj(TKappaPr)*KappaPr*TSigmax - 40*Conj(KappaPr)*((4*mphi2 + ms2
      + msbar2)*AbsSqr(Sigmax)*KappaPr + Conj(TSigmax)*(Sigmax*TKappaPr +
      KappaPr*TSigmax))));


   return beta_ms2;
}

} // namespace flexiblesusy
