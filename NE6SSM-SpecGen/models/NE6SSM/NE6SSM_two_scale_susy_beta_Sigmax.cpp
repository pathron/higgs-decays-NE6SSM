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

// File generated at Fri 20 Jun 2014 21:37:03

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Sigmax.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_Sigmax_one_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_Sigmax;

   beta_Sigmax = oneOver16PiSqr*(3*traceKappaAdjKappa*Sigmax + 2*AbsSqr(
      KappaPr)*Sigmax + 2*AbsSqr(Lambdax)*Sigmax - 0.1*Sigmax*Sqr(g1p)*Sqr(QS)
      + 3*Conj(Sigmax)*Sqr(Sigmax));


   return beta_Sigmax;
}

/**
 * Calculates the two-loop beta function of Sigmax.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_Sigmax_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;


   double beta_Sigmax;

   beta_Sigmax = -0.0025*twoLoop*Sigmax*(-4*Power(g1p,4)*Power(QS,4) +
      2400*traceKappaAdjKappaKappaAdjKappa + 2400*traceKappaAdjKappa*AbsSqr(
      Sigmax) + 3200*AbsSqr(KappaPr)*AbsSqr(Sigmax) - 320*traceKappaAdjKappa*
      Sqr(g1) - 780*traceKappaAdjKappa*Sqr(g1p) - 6400*traceKappaAdjKappa*Sqr(
      g3) - 249*Power(g1p,4)*Sqr(QS) + 60*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS) -
      40*AbsSqr(Sigmax)*Sqr(g1p)*Sqr(QS) - 40*AbsSqr(Lambdax)*(-60*
      traceYdAdjYd - 20*traceYeAdjYe - 60*traceYuAdjYu - 40*AbsSqr(Sigmax) + 12
      *Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2) - Sqr(g1p)*Sqr(QS)) + 3200*Sqr(Conj(
      KappaPr))*Sqr(KappaPr) + 1600*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 2400*Sqr(
      Conj(Sigmax))*Sqr(Sigmax));


   return beta_Sigmax;
}

} // namespace flexiblesusy
