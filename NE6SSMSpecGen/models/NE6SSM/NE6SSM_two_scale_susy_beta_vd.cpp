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

// File generated at Fri 13 Jun 2014 15:35:57

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vd.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = 0.05*oneOver16PiSqr*vd*(-60*traceYdAdjYd - 20*traceYeAdjYe -
      20*AbsSqr(Lambdax) + 6*Sqr(g1) + 9*Sqr(g1p) + 30*Sqr(g2));


   return beta_vd;
}

/**
 * Calculates the two-loop beta function of vd.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vd_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_vd;

   beta_vd = -0.000625*twoLoop*vd*(2088*Power(g1,4) + 2403*Power(g1p,4) +
      3400*Power(g2,4) - 14400*traceYdAdjYdYdAdjYd - 4800*traceYdAdjYuYuAdjYd
      - 4800*traceYeAdjYeYeAdjYe + 2400*traceYeAdjYe*Sqr(g1) + 400*traceYeAdjYe
      *Sqr(g1p) - 216*Sqr(g1)*Sqr(g1p) + 2400*traceYeAdjYe*Sqr(g2) + 720*Sqr(g1
      )*Sqr(g2) + 1080*Sqr(g1p)*Sqr(g2) + 400*traceYdAdjYd*(2*Sqr(g1) + 3*Sqr(
      g1p) + 18*Sqr(g2) + 64*Sqr(g3)) + 18*Power(g1p,4)*Sqr(QS) + 80*AbsSqr(
      Lambdax)*(-60*traceKappaAdjKappa - 60*traceYuAdjYu - 20*AbsSqr(Sigmax) +
      6*Sqr(g1) + 4*Sqr(g1p) + 30*Sqr(g2) + Sqr(g1p)*Sqr(QS)) - 4800*Sqr(Conj(
      Lambdax))*Sqr(Lambdax));


   return beta_vd;
}

} // namespace flexiblesusy
