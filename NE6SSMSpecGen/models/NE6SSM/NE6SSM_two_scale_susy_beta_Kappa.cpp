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
 * Calculates the one-loop beta function of Kappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_susy_parameters::calc_beta_Kappa_one_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = oneOver16PiSqr*(Kappa*(3*traceKappaAdjKappa + 2*AbsSqr(
      Lambdax) + AbsSqr(Sigmax) - 0.26666666666666666*Sqr(g1) - 0.65*Sqr(g1p) -
      5.333333333333333*Sqr(g3) - 0.05*Sqr(g1p)*Sqr(QS)) + 2*(Kappa*(Kappa)
      .adjoint()*Kappa));


   return beta_Kappa;
}

/**
 * Calculates the two-loop beta function of Kappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_susy_parameters::calc_beta_Kappa_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = twoLoop*(Kappa*(2.2755555555555556*Power(g1,4) + 4.28875*
      Power(g1p,4) + 14.222222222222221*Power(g3,4) + 0.005*Power(g1p,4)*Power(
      QS,4) - 6*traceKappaAdjKappaKappaAdjKappa - 2*AbsSqr(KappaPr)*AbsSqr(
      Sigmax) + 0.8*traceKappaAdjKappa*Sqr(g1) + 1.95*traceKappaAdjKappa*Sqr(
      g1p) + 0.29333333333333333*Sqr(g1)*Sqr(g1p) + 16*traceKappaAdjKappa*Sqr(
      g3) + 1.4222222222222223*Sqr(g1)*Sqr(g3) + 3.466666666666667*Sqr(g1p)*Sqr
      (g3) + 0.34375*Power(g1p,4)*Sqr(QS) - 0.15*traceKappaAdjKappa*Sqr(g1p)*
      Sqr(QS) + 0.1*AbsSqr(Lambdax)*(-60*traceYdAdjYd - 20*traceYeAdjYe - 60*
      traceYuAdjYu + 12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2) - Sqr(g1p)*Sqr(QS))
      - 4*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax)) +
      (-6*traceKappaAdjKappa - 4*AbsSqr(Lambdax) - 2*AbsSqr(Sigmax) + 0.1*Sqr(
      g1p)*Sqr(QS))*(Kappa*(Kappa).adjoint()*Kappa) - 2*(Kappa*(Kappa).adjoint(
      )*Kappa*(Kappa).adjoint()*Kappa));


   return beta_Kappa;
}

} // namespace flexiblesusy
