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

// File generated at Fri 20 Jun 2014 21:37:04

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vs.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vs_one_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_vs;

   beta_vs = oneOver16PiSqr*(-3*traceKappaAdjKappa*vs - 2*vs*AbsSqr(
      Lambdax) - vs*AbsSqr(Sigmax) + 0.05*vs*Sqr(g1p)*Sqr(QS));


   return beta_vs;
}

/**
 * Calculates the two-loop beta function of vs.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vs_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;


   double beta_vs;

   beta_vs = -0.000625*twoLoop*vs*(4*Power(g1p,4)*Power(QS,4) - 9600*
      traceKappaAdjKappaKappaAdjKappa + 1280*traceKappaAdjKappa*Sqr(g1) + 3120*
      traceKappaAdjKappa*Sqr(g1p) + 160*AbsSqr(Lambdax)*(-60*traceYdAdjYd - 20*
      traceYeAdjYe - 60*traceYuAdjYu + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2)) +
      25600*traceKappaAdjKappa*Sqr(g3) + 249*Power(g1p,4)*Sqr(QS) + 80*AbsSqr(
      Sigmax)*(-40*AbsSqr(KappaPr) + Sqr(g1p)*Sqr(QS)) - 6400*Sqr(Conj(Lambdax)
      )*Sqr(Lambdax) - 3200*Sqr(Conj(Sigmax))*Sqr(Sigmax));


   return beta_vs;
}

} // namespace flexiblesusy
