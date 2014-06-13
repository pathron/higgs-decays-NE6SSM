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

// File generated at Fri 13 Jun 2014 15:35:55

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 0.4666666666666667*Sqr(g1) - 0.7*Sqr(g1p) - 3*Sqr(g2) -
      5.333333333333333*Sqr(g3)) + 3*(Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu)
      ;


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = twoLoop*(Yd*(4.028888888888889*Power(g1,4) + 4.6025*Power(
      g1p,4) + 10.5*Power(g2,4) + 14.222222222222221*Power(g3,4) - 9*
      traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe + 1.2
      *traceYeAdjYe*Sqr(g1) - 0.2*traceYeAdjYe*Sqr(g1p) - 0.5133333333333333*
      Sqr(g1)*Sqr(g1p) + Sqr(g1)*Sqr(g2) + 1.5*Sqr(g1p)*Sqr(g2) - 0.2*
      traceYdAdjYd*(2*Sqr(g1) + 3*Sqr(g1p) - 80*Sqr(g3)) + 0.8888888888888888*
      Sqr(g1)*Sqr(g3) + 1.3333333333333333*Sqr(g1p)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3)
      + 0.035*Power(g1p,4)*Sqr(QS) + 0.05*AbsSqr(Lambdax)*(-60*
      traceKappaAdjKappa - 60*traceYuAdjYu - 20*AbsSqr(Sigmax) - 5*Sqr(g1p) +
      Sqr(g1p)*Sqr(QS)) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd
      - 3*traceYeAdjYe - 3*AbsSqr(Lambdax) + 0.8*Sqr(g1) + 1.2*Sqr(g1p) + 6*
      Sqr(g2))*(Yd*Yd.adjoint()*Yd) - 3*traceYuAdjYu*(Yd*Yu.adjoint()*Yu) -
      AbsSqr(Lambdax)*(Yd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(Yd*Yu.adjoint()*Yu) +
      0.2*Sqr(g1p)*(Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*
      Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu));


   return beta_Yd;
}

} // namespace flexiblesusy
