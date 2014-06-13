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
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = 0.1*oneOver16PiSqr*vu*(-30*traceYuAdjYu - 10*AbsSqr(Lambdax)
      + 3*Sqr(g1) + 2*Sqr(g1p) + 15*Sqr(g2));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = -0.0025*twoLoop*vu*(522*Power(g1,4) + 257*Power(g1p,4) + 850
      *Power(g2,4) - 1200*traceYdAdjYuYuAdjYd - 3600*traceYuAdjYuYuAdjYu + 96*
      Sqr(g1)*Sqr(g1p) + 180*Sqr(g1)*Sqr(g2) + 120*Sqr(g1p)*Sqr(g2) + 40*
      traceYuAdjYu*(17*Sqr(g1) + 3*Sqr(g1p) + 45*Sqr(g2) + 160*Sqr(g3)) + 2*
      Power(g1p,4)*Sqr(QS) + 20*AbsSqr(Lambdax)*(-60*traceKappaAdjKappa - 60*
      traceYdAdjYd - 20*traceYeAdjYe - 20*AbsSqr(Sigmax) + 6*Sqr(g1) + 9*Sqr(
      g1p) + 30*Sqr(g2) + Sqr(g1p)*Sqr(QS)) - 1200*Sqr(Conj(Lambdax))*Sqr(
      Lambdax));


   return beta_vu;
}

} // namespace flexiblesusy
