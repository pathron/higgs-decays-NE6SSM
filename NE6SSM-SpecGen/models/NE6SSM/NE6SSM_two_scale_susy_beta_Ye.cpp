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

// File generated at Fri 20 Jun 2014 21:37:02

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 0.7*Sqr(g1p) - 3*Sqr(g2)) + 3*(Ye*Ye.adjoint()*
      Ye));


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = twoLoop*(0.0025*Ye*(6696*Power(g1,4) + 1841*Power(g1p,4) +
      4200*Power(g2,4) - 3600*traceYdAdjYdYdAdjYd - 1200*traceYdAdjYuYuAdjYd -
      1200*traceYeAdjYeYeAdjYe + 480*traceYeAdjYe*Sqr(g1) - 80*traceYeAdjYe*Sqr
      (g1p) - 12*Sqr(g1)*Sqr(g1p) + 720*Sqr(g1)*Sqr(g2) + 780*Sqr(g1p)*Sqr(g2)
      - 80*traceYdAdjYd*(2*Sqr(g1) + 3*Sqr(g1p) - 80*Sqr(g3)) + 14*Power(g1p,4)
      *Sqr(QS) + 20*AbsSqr(Lambdax)*(-60*traceKappaAdjKappa - 60*traceYuAdjYu -
      20*AbsSqr(Sigmax) - 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS)) - 1200*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)) + 1.5*(-6*traceYdAdjYd - 2*traceYeAdjYe - 2*
      AbsSqr(Lambdax) + Sqr(g1p) + 4*Sqr(g2))*(Ye*Ye.adjoint()*Ye) - 4*(Ye*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye));


   return beta_Ye;
}

} // namespace flexiblesusy
