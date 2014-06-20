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
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_Lambdax;

   beta_Lambdax = 0.05*oneOver16PiSqr*Lambdax*(60*traceKappaAdjKappa + 60
      *traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu + 80*AbsSqr(Lambdax) +
      20*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2) - Sqr(g1p)*Sqr(
      QS));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;


   double beta_Lambdax;

   beta_Lambdax = -0.00125*twoLoop*Lambdax*(-4176*Power(g1,4) - 3431*
      Power(g1p,4) - 8400*Power(g2,4) - 4*Power(g1p,4)*Power(QS,4) + 4800*
      traceKappaAdjKappaKappaAdjKappa + 7200*traceYdAdjYdYdAdjYd + 4800*
      traceYdAdjYuYuAdjYd + 2400*traceYeAdjYeYeAdjYe + 7200*traceYuAdjYuYuAdjYu
      + 1600*AbsSqr(KappaPr)*AbsSqr(Sigmax) - 640*traceKappaAdjKappa*Sqr(g1) +
      320*traceYdAdjYd*Sqr(g1) - 960*traceYeAdjYe*Sqr(g1) - 640*traceYuAdjYu*
      Sqr(g1) - 1560*traceKappaAdjKappa*Sqr(g1p) + 480*traceYdAdjYd*Sqr(g1p) +
      160*traceYeAdjYe*Sqr(g1p) + 240*traceYuAdjYu*Sqr(g1p) - 168*Sqr(g1)*Sqr(
      g1p) - 1440*Sqr(g1)*Sqr(g2) - 1560*Sqr(g1p)*Sqr(g2) - 80*AbsSqr(Lambdax)*
      (-60*traceKappaAdjKappa - 90*traceYdAdjYd - 30*traceYeAdjYe - 90*
      traceYuAdjYu - 20*AbsSqr(Sigmax) + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2))
      - 12800*traceKappaAdjKappa*Sqr(g3) - 12800*traceYdAdjYd*Sqr(g3) - 12800*
      traceYuAdjYu*Sqr(g3) - 275*Power(g1p,4)*Sqr(QS) + 120*traceKappaAdjKappa*
      Sqr(g1p)*Sqr(QS) + 8000*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 1600*Sqr(Conj(
      Sigmax))*Sqr(Sigmax));


   return beta_Lambdax;
}

} // namespace flexiblesusy
