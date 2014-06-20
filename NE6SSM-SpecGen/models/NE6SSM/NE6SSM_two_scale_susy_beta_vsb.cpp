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
 * Calculates the one-loop beta function of vsb.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vsb_one_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);


   double beta_vsb;

   beta_vsb = oneOver16PiSqr*(-(vsb*AbsSqr(Sigmax)) + 0.05*vsb*Sqr(g1p)*
      Sqr(QS));


   return beta_vsb;
}

/**
 * Calculates the two-loop beta function of vsb.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vsb_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_vsb;

   beta_vsb = -0.000625*twoLoop*vsb*(Power(g1p,4)*Sqr(QS)*(249 + 4*Sqr(QS
      )) + 80*AbsSqr(Sigmax)*(-60*traceKappaAdjKappa - 40*AbsSqr(KappaPr) - 40*
      AbsSqr(Lambdax) + Sqr(g1p)*Sqr(QS)) - 3200*Sqr(Conj(Sigmax))*Sqr(Sigmax))
      ;


   return beta_vsb;
}

} // namespace flexiblesusy
