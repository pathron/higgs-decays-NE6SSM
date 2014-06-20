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

// File generated at Fri 20 Jun 2014 21:37:49

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 -
      0.9486832980505138*g1p*Tr14 + 6*traceconjTYdTpTYd + 2*traceconjTYeTpTYe +
      6*tracemd2YdAdjYd + 2*traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*
      tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*mHd2*
      AbsSqr(Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*
      AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 1.8*AbsSqr(MassBp)*Sqr(g1p
      ) - 6*AbsSqr(MassWB)*Sqr(g2));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   beta_mHd2 = 0.005*twoLoop*(4*Conj(MassB)*Sqr(g1)*(40*traceAdjYdTYd -
      120*traceAdjYeTYe - 80*MassB*traceYdAdjYd + 240*MassB*traceYeAdjYe + 1566
      *MassB*Sqr(g1) - 54*MassB*Sqr(g1p) - 27*MassBp*Sqr(g1p) + 180*MassB*Sqr(
      g2) + 90*MassWB*Sqr(g2)) + Conj(MassBp)*Sqr(g1p)*(240*traceAdjYdTYd + 80*
      traceAdjYeTYe - 480*MassBp*traceYdAdjYd - 160*MassBp*traceYeAdjYe - 108*
      MassB*Sqr(g1) - 216*MassBp*Sqr(g1) + 7209*MassBp*Sqr(g1p) + 1080*MassBp*
      Sqr(g2) + 540*MassWB*Sqr(g2) + 54*MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(
      Lambdax)*(-5 + Sqr(QS))*(2*MassBp*Lambdax - TLambdax)) + 20*(60*Power(g2,
      4)*Tr22 + 14.696938456699067*g1*g1p*Tr2U114 + 14.696938456699067*g1*g1p*
      Tr2U141 - 30.983866769659336*g1*Tr31 - 37.94733192202055*g1p*Tr34 - 360*
      tracemd2YdAdjYdYdAdjYd - 60*tracemd2YdAdjYuYuAdjYd - 120*
      traceme2YeAdjYeYeAdjYe - 120*traceml2AdjYeYeAdjYeYe - 360*
      tracemq2AdjYdYdAdjYdYd - 60*tracemq2AdjYdYdAdjYuYu - 60*
      tracemq2AdjYuYuAdjYdYd - 60*tracemu2YuAdjYdYdAdjYu - 360*
      traceYdAdjTYdTYdAdjYd - 60*traceYdAdjTYuTYuAdjYd - 360*
      traceYdAdjYdTYdAdjTYd - 360*mHd2*traceYdAdjYdYdAdjYd - 60*
      traceYdAdjYuTYuAdjTYd - 60*mHd2*traceYdAdjYuYuAdjYd - 60*mHu2*
      traceYdAdjYuYuAdjYd - 120*traceYeAdjTYeTYeAdjYe - 120*
      traceYeAdjYeTYeAdjTYe - 120*mHd2*traceYeAdjYeYeAdjYe - 60*
      traceYuAdjTYdTYdAdjYu - 60*traceYuAdjYdTYdAdjTYu - 60*traceKappaAdjKappa*
      AbsSqr(TLambdax) - 60*traceYuAdjYu*AbsSqr(TLambdax) - 20*AbsSqr(Sigmax)*
      AbsSqr(TLambdax) - 60*traceAdjKappaTKappa*Conj(TLambdax)*Lambdax - 60*
      traceAdjYuTYu*Conj(TLambdax)*Lambdax + 12*Tr2U111*Sqr(g1) - 8*
      traceconjTYdTpTYd*Sqr(g1) + 8*MassB*traceconjTYdTpYd*Sqr(g1) + 24*
      traceconjTYeTpTYe*Sqr(g1) - 24*MassB*traceconjTYeTpYe*Sqr(g1) - 8*
      tracemd2YdAdjYd*Sqr(g1) + 24*traceme2YeAdjYe*Sqr(g1) + 24*traceml2AdjYeYe
      *Sqr(g1) - 8*tracemq2AdjYdYd*Sqr(g1) - 8*mHd2*traceYdAdjYd*Sqr(g1) + 24*
      mHd2*traceYeAdjYe*Sqr(g1) + 18*Tr2U144*Sqr(g1p) - 12*traceconjTYdTpTYd*
      Sqr(g1p) + 12*MassBp*traceconjTYdTpYd*Sqr(g1p) - 4*traceconjTYeTpTYe*Sqr(
      g1p) + 4*MassBp*traceconjTYeTpYe*Sqr(g1p) - 12*tracemd2YdAdjYd*Sqr(g1p) -
      4*traceme2YeAdjYe*Sqr(g1p) - 4*traceml2AdjYeYe*Sqr(g1p) - 12*
      tracemq2AdjYdYd*Sqr(g1p) - 12*mHd2*traceYdAdjYd*Sqr(g1p) - 4*mHd2*
      traceYeAdjYe*Sqr(g1p) - 5*AbsSqr(TLambdax)*Sqr(g1p) + 5*MassBp*Conj(
      TLambdax)*Lambdax*Sqr(g1p) + 3*Conj(MassWB)*Sqr(g2)*(6*(MassB + 2*MassWB)
      *Sqr(g1) + 9*(MassBp + 2*MassWB)*Sqr(g1p) + 170*MassWB*Sqr(g2)) + 320*
      traceconjTYdTpTYd*Sqr(g3) - 320*MassG*traceconjTYdTpYd*Sqr(g3) + 320*
      tracemd2YdAdjYd*Sqr(g3) + 320*tracemq2AdjYdYd*Sqr(g3) + 320*mHd2*
      traceYdAdjYd*Sqr(g3) + 640*traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 320*
      traceAdjYdTYd*Conj(MassG)*Sqr(g3) + AbsSqr(TLambdax)*Sqr(g1p)*Sqr(QS) -
      MassBp*Conj(TLambdax)*Lambdax*Sqr(g1p)*Sqr(QS) - 120*(mHd2 + mHu2 + ms2)*
      Sqr(Conj(Lambdax))*Sqr(Lambdax) + Conj(Lambdax)*(-60*
      traceconjTKappaTpTKappa*Lambdax - 60*traceconjTYuTpTYu*Lambdax - 60*mHd2*
      traceKappaAdjKappa*Lambdax - 60*mHu2*traceKappaAdjKappa*Lambdax - 120*ms2
      *traceKappaAdjKappa*Lambdax - 60*traceKappaAdjKappaconjmDx2*Lambdax - 60*
      traceKappaconjmDxbar2AdjKappa*Lambdax - 60*tracemq2AdjYuYu*Lambdax - 60*
      tracemu2YuAdjYu*Lambdax - 60*mHd2*traceYuAdjYu*Lambdax - 120*mHu2*
      traceYuAdjYu*Lambdax - 60*ms2*traceYuAdjYu*Lambdax - 20*(mHd2 + mHu2 +
      mphi2 + 2*ms2 + msbar2)*AbsSqr(Sigmax)*Lambdax - 240*AbsSqr(TLambdax)*
      Lambdax - 20*AbsSqr(TSigmax)*Lambdax - 5*mHd2*Lambdax*Sqr(g1p) - 5*mHu2*
      Lambdax*Sqr(g1p) - 5*ms2*Lambdax*Sqr(g1p) + mHd2*Lambdax*Sqr(g1p)*Sqr(QS)
      + mHu2*Lambdax*Sqr(g1p)*Sqr(QS) + ms2*Lambdax*Sqr(g1p)*Sqr(QS) - 60*
      traceconjTKappaTpKappa*TLambdax - 60*traceconjTYuTpYu*TLambdax - 20*Conj(
      TSigmax)*Sigmax*TLambdax) - 20*Conj(Sigmax)*Conj(TLambdax)*Lambdax*
      TSigmax));


   return beta_mHd2;
}

} // namespace flexiblesusy
