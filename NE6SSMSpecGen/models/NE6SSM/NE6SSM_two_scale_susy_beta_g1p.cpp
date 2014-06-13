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
 * Calculates the one-loop beta function of g1p.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_g1p_one_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);


   double beta_g1p;

   beta_g1p = 0.025*Power(g1p,3)*oneOver16PiSqr*(249 + 2*Sqr(QS));


   return beta_g1p;
}

/**
 * Calculates the two-loop beta function of g1p.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_g1p_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_g1p;

   beta_g1p = 0.0025*Power(g1p,3)*twoLoop*(-780*traceKappaAdjKappa - 1680
      *traceYdAdjYd - 560*traceYeAdjYe - 720*traceYuAdjYu + 984*Sqr(g1) + 1401*
      Sqr(g1p) + 2*Power(QS,4)*Sqr(g1p) + 2520*Sqr(g2) + 9600*Sqr(g3) - 60*
      traceKappaAdjKappa*Sqr(QS) - 40*AbsSqr(Sigmax)*Sqr(QS) - 40*AbsSqr(
      Lambdax)*(13 + Sqr(QS)));


   return beta_g1p;
}

} // namespace flexiblesusy
