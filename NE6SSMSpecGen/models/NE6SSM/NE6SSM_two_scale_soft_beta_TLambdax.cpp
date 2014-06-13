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

// File generated at Fri 13 Jun 2014 15:36:38

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_TLambdax;

   beta_TLambdax = oneOver16PiSqr*((3*traceKappaAdjKappa + 3*traceYdAdjYd
      + traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) + AbsSqr(Sigmax) -
      0.6*Sqr(g1) - 0.65*Sqr(g1p) - 3*Sqr(g2) - 0.05*Sqr(g1p)*Sqr(QS))*TLambdax
      + 0.1*Lambdax*(60*traceAdjKappaTKappa + 60*traceAdjYdTYd + 20*
      traceAdjYeTYe + 60*traceAdjYuTYu + 12*MassB*Sqr(g1) + 13*MassBp*Sqr(g1p)
      + 60*MassWB*Sqr(g2) + MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))
      ;


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;


   double beta_TLambdax;

   beta_TLambdax = twoLoop*(-20.88*Power(g1,4)*MassB*Lambdax - 17.155*
      Power(g1p,4)*MassBp*Lambdax - 42*Power(g2,4)*MassWB*Lambdax - 0.02*Power(
      g1p,4)*MassBp*Power(QS,4)*Lambdax - 24*traceKappaAdjKappaTKappaAdjKappa*
      Lambdax - 36*traceYdAdjYdTYdAdjYd*Lambdax - 12*traceYdAdjYuTYuAdjYd*
      Lambdax - 12*traceYeAdjYeTYeAdjYe*Lambdax - 12*traceYuAdjYdTYdAdjYu*
      Lambdax - 36*traceYuAdjYuTYuAdjYu*Lambdax + 1.6*traceAdjKappaTKappa*
      Lambdax*Sqr(g1) - 0.8*traceAdjYdTYd*Lambdax*Sqr(g1) + 2.4*traceAdjYeTYe*
      Lambdax*Sqr(g1) + 1.6*traceAdjYuTYu*Lambdax*Sqr(g1) - 1.6*MassB*
      traceKappaAdjKappa*Lambdax*Sqr(g1) + 0.8*MassB*traceYdAdjYd*Lambdax*Sqr(
      g1) - 2.4*MassB*traceYeAdjYe*Lambdax*Sqr(g1) - 1.6*MassB*traceYuAdjYu*
      Lambdax*Sqr(g1) + 3.9*traceAdjKappaTKappa*Lambdax*Sqr(g1p) - 1.2*
      traceAdjYdTYd*Lambdax*Sqr(g1p) - 0.4*traceAdjYeTYe*Lambdax*Sqr(g1p) - 0.6
      *traceAdjYuTYu*Lambdax*Sqr(g1p) - 3.9*MassBp*traceKappaAdjKappa*Lambdax*
      Sqr(g1p) + 1.2*MassBp*traceYdAdjYd*Lambdax*Sqr(g1p) + 0.4*MassBp*
      traceYeAdjYe*Lambdax*Sqr(g1p) + 0.6*MassBp*traceYuAdjYu*Lambdax*Sqr(g1p)
      - 0.78*MassB*Lambdax*Sqr(g1)*Sqr(g1p) - 0.78*MassBp*Lambdax*Sqr(g1)*Sqr(
      g1p) - 3.6*MassB*Lambdax*Sqr(g1)*Sqr(g2) - 3.6*MassWB*Lambdax*Sqr(g1)*Sqr
      (g2) - 3.9*MassBp*Lambdax*Sqr(g1p)*Sqr(g2) - 3.9*MassWB*Lambdax*Sqr(g1p)*
      Sqr(g2) + 32*traceAdjKappaTKappa*Lambdax*Sqr(g3) + 32*traceAdjYdTYd*
      Lambdax*Sqr(g3) + 32*traceAdjYuTYu*Lambdax*Sqr(g3) - 32*MassG*
      traceKappaAdjKappa*Lambdax*Sqr(g3) - 32*MassG*traceYdAdjYd*Lambdax*Sqr(g3
      ) - 32*MassG*traceYuAdjYu*Lambdax*Sqr(g3) - 1.375*Power(g1p,4)*MassBp*
      Lambdax*Sqr(QS) - 0.3*traceAdjKappaTKappa*Lambdax*Sqr(g1p)*Sqr(QS) + 0.3*
      MassBp*traceKappaAdjKappa*Lambdax*Sqr(g1p)*Sqr(QS) + 5.22*Power(g1,4)*
      TLambdax + 4.28875*Power(g1p,4)*TLambdax + 10.5*Power(g2,4)*TLambdax +
      0.005*Power(g1p,4)*Power(QS,4)*TLambdax - 6*
      traceKappaAdjKappaKappaAdjKappa*TLambdax - 9*traceYdAdjYdYdAdjYd*TLambdax
      - 6*traceYdAdjYuYuAdjYd*TLambdax - 3*traceYeAdjYeYeAdjYe*TLambdax - 9*
      traceYuAdjYuYuAdjYu*TLambdax + 0.8*traceKappaAdjKappa*Sqr(g1)*TLambdax -
      0.4*traceYdAdjYd*Sqr(g1)*TLambdax + 1.2*traceYeAdjYe*Sqr(g1)*TLambdax +
      0.8*traceYuAdjYu*Sqr(g1)*TLambdax + 1.95*traceKappaAdjKappa*Sqr(g1p)*
      TLambdax - 0.6*traceYdAdjYd*Sqr(g1p)*TLambdax - 0.2*traceYeAdjYe*Sqr(g1p)
      *TLambdax - 0.3*traceYuAdjYu*Sqr(g1p)*TLambdax + 0.39*Sqr(g1)*Sqr(g1p)*
      TLambdax + 1.8*Sqr(g1)*Sqr(g2)*TLambdax + 1.95*Sqr(g1p)*Sqr(g2)*TLambdax
      + 16*traceKappaAdjKappa*Sqr(g3)*TLambdax + 16*traceYdAdjYd*Sqr(g3)*
      TLambdax + 16*traceYuAdjYu*Sqr(g3)*TLambdax + 0.34375*Power(g1p,4)*Sqr(QS
      )*TLambdax - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS)*TLambdax - 50*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)*TLambdax - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)*
      TLambdax - 8*Lambdax*Sigmax*Sqr(Conj(Sigmax))*TSigmax - 2*Conj(KappaPr)*
      Conj(Sigmax)*(2*Lambdax*Sigmax*TKappaPr + KappaPr*Sigmax*TLambdax + 2*
      KappaPr*Lambdax*TSigmax) - 0.1*AbsSqr(Lambdax)*(-3*(-60*
      traceKappaAdjKappa - 90*traceYdAdjYd - 30*traceYeAdjYe - 90*traceYuAdjYu
      - 20*AbsSqr(Sigmax) + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2))*TLambdax + 2
      *Lambdax*(60*traceAdjKappaTKappa + 90*traceAdjYdTYd + 30*traceAdjYeTYe +
      90*traceAdjYuTYu + 12*MassB*Sqr(g1) + 13*MassBp*Sqr(g1p) + 60*MassWB*Sqr(
      g2) + 20*Conj(Sigmax)*TSigmax)));


   return beta_TLambdax;
}

} // namespace flexiblesusy
