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

// File generated at Fri 13 Jun 2014 15:35:56

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g1.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_g1_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = 8.4*Power(g1,3)*oneOver16PiSqr;


   return beta_g1;
}

/**
 * Calculates the two-loop beta function of g1.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_g1_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_g1;

   beta_g1 = 0.02*Power(g1,3)*twoLoop*(-40*traceKappaAdjKappa - 140*
      traceYdAdjYd - 180*traceYeAdjYe - 260*traceYuAdjYu - 60*AbsSqr(Lambdax) +
      432*Sqr(g1) + 123*Sqr(g1p) + 360*Sqr(g2) + 1200*Sqr(g3));


   return beta_g1;
}

} // namespace flexiblesusy
