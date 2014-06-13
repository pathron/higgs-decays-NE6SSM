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

// File generated at Fri 13 Jun 2014 15:36:36

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = oneOver16PiSqr*(3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd - 0.7*Sqr(g1p)*TYd -
      3*Sqr(g2)*TYd - 5.333333333333333*Sqr(g3)*TYd + Yd*(6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 0.9333333333333333*MassB*Sqr(g1) + 1.4*MassBp*Sqr(g1p) +
      6*MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 2*Conj(Lambdax)*
      TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*
      Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu);


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = twoLoop*(4.028888888888889*Power(g1,4)*TYd + 4.6025*Power(
      g1p,4)*TYd + 10.5*Power(g2,4)*TYd + 14.222222222222221*Power(g3,4)*TYd -
      9*traceYdAdjYdYdAdjYd*TYd - 3*traceYdAdjYuYuAdjYd*TYd - 3*
      traceYeAdjYeYeAdjYe*TYd - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*TYd - 3*
      traceYuAdjYu*AbsSqr(Lambdax)*TYd - AbsSqr(Lambdax)*AbsSqr(Sigmax)*TYd -
      0.4*traceYdAdjYd*Sqr(g1)*TYd + 1.2*traceYeAdjYe*Sqr(g1)*TYd - 0.6*
      traceYdAdjYd*Sqr(g1p)*TYd - 0.2*traceYeAdjYe*Sqr(g1p)*TYd - 0.25*AbsSqr(
      Lambdax)*Sqr(g1p)*TYd + 0.32666666666666666*Sqr(g1)*Sqr(g1p)*TYd + Sqr(g1
      )*Sqr(g2)*TYd + 1.5*Sqr(g1p)*Sqr(g2)*TYd + 16*traceYdAdjYd*Sqr(g3)*TYd +
      0.8888888888888888*Sqr(g1)*Sqr(g3)*TYd + 1.3333333333333333*Sqr(g1p)*Sqr(
      g3)*TYd + 8*Sqr(g2)*Sqr(g3)*TYd + 0.035*Power(g1p,4)*Sqr(QS)*TYd + 0.05*
      AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TYd - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*
      TYd - 0.0011111111111111111*Yd*(14504*Power(g1,4)*MassB + 16569*Power(g1p
      ,4)*MassBp + 51200*Power(g3,4)*MassG + 37800*Power(g2,4)*MassWB + 32400*
      traceYdAdjYdTYdAdjYd + 5400*traceYdAdjYuTYuAdjYd + 10800*
      traceYeAdjYeTYeAdjYe + 5400*traceYuAdjYdTYdAdjYu + 720*traceAdjYdTYd*Sqr(
      g1) - 2160*traceAdjYeTYe*Sqr(g1) + 2160*MassB*traceYeAdjYe*Sqr(g1) + 1080
      *traceAdjYdTYd*Sqr(g1p) + 360*traceAdjYeTYe*Sqr(g1p) - 360*MassBp*
      traceYeAdjYe*Sqr(g1p) + 588*MassB*Sqr(g1)*Sqr(g1p) + 588*MassBp*Sqr(g1)*
      Sqr(g1p) + 1800*MassB*Sqr(g1)*Sqr(g2) + 1800*MassWB*Sqr(g1)*Sqr(g2) +
      2700*MassBp*Sqr(g1p)*Sqr(g2) + 2700*MassWB*Sqr(g1p)*Sqr(g2) - 28800*
      traceAdjYdTYd*Sqr(g3) + 1600*MassB*Sqr(g1)*Sqr(g3) + 1600*MassG*Sqr(g1)*
      Sqr(g3) + 2400*MassBp*Sqr(g1p)*Sqr(g3) + 2400*MassG*Sqr(g1p)*Sqr(g3) +
      14400*MassG*Sqr(g2)*Sqr(g3) + 14400*MassWB*Sqr(g2)*Sqr(g3) - 360*
      traceYdAdjYd*(2*MassB*Sqr(g1) + 3*MassBp*Sqr(g1p) - 80*MassG*Sqr(g3)) +
      126*Power(g1p,4)*MassBp*Sqr(QS) + 10800*Lambdax*Sqr(Conj(Lambdax))*
      TLambdax + 90*Conj(Lambdax)*((60*traceKappaAdjKappa + 60*traceYuAdjYu +
      20*AbsSqr(Sigmax) + 5*Sqr(g1p) - Sqr(g1p)*Sqr(QS))*TLambdax + Lambdax*(60
      *traceAdjKappaTKappa + 60*traceAdjYuTYu - 5*MassBp*Sqr(g1p) + MassBp*Sqr(
      g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) - 0.4*(45*traceAdjYdTYd + 15*
      traceAdjYeTYe + 4*MassB*Sqr(g1) + 6*MassBp*Sqr(g1p) + 30*MassWB*Sqr(g2) +
      15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) - 12*traceYdAdjYd*(Yd*
      Yd.adjoint()*TYd) - 4*traceYeAdjYe*(Yd*Yd.adjoint()*TYd) - 4*AbsSqr(
      Lambdax)*(Yd*Yd.adjoint()*TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint()*TYd) + 1.8*
      Sqr(g1p)*(Yd*Yd.adjoint()*TYd) + 6*Sqr(g2)*(Yd*Yd.adjoint()*TYd) - 6*
      traceAdjYuTYu*(Yd*Yu.adjoint()*Yu) - 1.6*MassB*Sqr(g1)*(Yd*Yu.adjoint()*
      Yu) - 0.4*MassBp*Sqr(g1p)*(Yd*Yu.adjoint()*Yu) - 2*Conj(Lambdax)*TLambdax
      *(Yd*Yu.adjoint()*Yu) - 6*traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 2*AbsSqr(
      Lambdax)*(Yd*Yu.adjoint()*TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu) + 0.4*
      Sqr(g1p)*(Yd*Yu.adjoint()*TYu) - 15*traceYdAdjYd*(TYd*Yd.adjoint()*Yd) -
      5*traceYeAdjYe*(TYd*Yd.adjoint()*Yd) - 5*AbsSqr(Lambdax)*(TYd*Yd.adjoint(
      )*Yd) + 1.2*Sqr(g1)*(TYd*Yd.adjoint()*Yd) + 1.8*Sqr(g1p)*(TYd*Yd.adjoint(
      )*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*Yd) - 3*traceYuAdjYu*(TYd*Yu.adjoint
      ()*Yu) - AbsSqr(Lambdax)*(TYd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(TYd*
      Yu.adjoint()*Yu) + 0.2*Sqr(g1p)*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint(
      )*Yd*Yd.adjoint()*TYd) - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*
      TYu) - 4*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*
      Yu.adjoint()*Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*
      Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu
      ));


   return beta_TYd;
}

} // namespace flexiblesusy
