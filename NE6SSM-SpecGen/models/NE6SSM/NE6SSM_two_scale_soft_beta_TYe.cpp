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

// File generated at Fri 20 Jun 2014 21:37:44

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = oneOver16PiSqr*(3*traceYdAdjYd*TYe + traceYeAdjYe*TYe +
      AbsSqr(Lambdax)*TYe - 1.8*Sqr(g1)*TYe - 0.7*Sqr(g1p)*TYe - 3*Sqr(g2)*TYe
      + Ye*(6*traceAdjYdTYd + 2*traceAdjYeTYe + 3.6*MassB*Sqr(g1) + 1.4*MassBp*
      Sqr(g1p) + 6*MassWB*Sqr(g2) + 2*Conj(Lambdax)*TLambdax) + 4*(Ye*
      Ye.adjoint()*TYe) + 5*(TYe*Ye.adjoint()*Ye));


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NE6SSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = twoLoop*(16.74*Power(g1,4)*TYe + 4.6025*Power(g1p,4)*TYe +
      10.5*Power(g2,4)*TYe - 9*traceYdAdjYdYdAdjYd*TYe - 3*traceYdAdjYuYuAdjYd*
      TYe - 3*traceYeAdjYeYeAdjYe*TYe - 3*traceKappaAdjKappa*AbsSqr(Lambdax)*
      TYe - 3*traceYuAdjYu*AbsSqr(Lambdax)*TYe - AbsSqr(Lambdax)*AbsSqr(Sigmax)
      *TYe - 0.4*traceYdAdjYd*Sqr(g1)*TYe + 1.2*traceYeAdjYe*Sqr(g1)*TYe - 0.6*
      traceYdAdjYd*Sqr(g1p)*TYe - 0.2*traceYeAdjYe*Sqr(g1p)*TYe - 0.25*AbsSqr(
      Lambdax)*Sqr(g1p)*TYe + 0.51*Sqr(g1)*Sqr(g1p)*TYe + 1.8*Sqr(g1)*Sqr(g2)*
      TYe + 1.95*Sqr(g1p)*Sqr(g2)*TYe + 16*traceYdAdjYd*Sqr(g3)*TYe + 0.035*
      Power(g1p,4)*Sqr(QS)*TYe + 0.05*AbsSqr(Lambdax)*Sqr(g1p)*Sqr(QS)*TYe - 3*
      Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYe - 0.01*Ye*(6696*Power(g1,4)*MassB +
      1841*Power(g1p,4)*MassBp + 4200*Power(g2,4)*MassWB + 3600*
      traceYdAdjYdTYdAdjYd + 600*traceYdAdjYuTYuAdjYd + 1200*
      traceYeAdjYeTYeAdjYe + 600*traceYuAdjYdTYdAdjYu + 80*traceAdjYdTYd*Sqr(g1
      ) - 240*traceAdjYeTYe*Sqr(g1) + 240*MassB*traceYeAdjYe*Sqr(g1) + 120*
      traceAdjYdTYd*Sqr(g1p) + 40*traceAdjYeTYe*Sqr(g1p) - 40*MassBp*
      traceYeAdjYe*Sqr(g1p) + 102*MassB*Sqr(g1)*Sqr(g1p) + 102*MassBp*Sqr(g1)*
      Sqr(g1p) + 360*MassB*Sqr(g1)*Sqr(g2) + 360*MassWB*Sqr(g1)*Sqr(g2) + 390*
      MassBp*Sqr(g1p)*Sqr(g2) + 390*MassWB*Sqr(g1p)*Sqr(g2) - 3200*
      traceAdjYdTYd*Sqr(g3) - 40*traceYdAdjYd*(2*MassB*Sqr(g1) + 3*MassBp*Sqr(
      g1p) - 80*MassG*Sqr(g3)) + 14*Power(g1p,4)*MassBp*Sqr(QS) + 1200*Lambdax*
      Sqr(Conj(Lambdax))*TLambdax + 10*Conj(Lambdax)*((60*traceKappaAdjKappa +
      60*traceYuAdjYu + 20*AbsSqr(Sigmax) + 5*Sqr(g1p) - Sqr(g1p)*Sqr(QS))*
      TLambdax + Lambdax*(60*traceAdjKappaTKappa + 60*traceAdjYuTYu - 5*MassBp*
      Sqr(g1p) + MassBp*Sqr(g1p)*Sqr(QS) + 20*Conj(Sigmax)*TSigmax))) - 3*(6*
      traceAdjYdTYd + 2*traceAdjYeTYe + MassBp*Sqr(g1p) + 4*MassWB*Sqr(g2) + 2*
      Conj(Lambdax)*TLambdax)*(Ye*Ye.adjoint()*Ye) - 12*traceYdAdjYd*(Ye*
      Ye.adjoint()*TYe) - 4*traceYeAdjYe*(Ye*Ye.adjoint()*TYe) - 4*AbsSqr(
      Lambdax)*(Ye*Ye.adjoint()*TYe) + 1.2*Sqr(g1)*(Ye*Ye.adjoint()*TYe) + 1.8*
      Sqr(g1p)*(Ye*Ye.adjoint()*TYe) + 6*Sqr(g2)*(Ye*Ye.adjoint()*TYe) - 15*
      traceYdAdjYd*(TYe*Ye.adjoint()*Ye) - 5*traceYeAdjYe*(TYe*Ye.adjoint()*Ye)
      - 5*AbsSqr(Lambdax)*(TYe*Ye.adjoint()*Ye) - 1.2*Sqr(g1)*(TYe*Ye.adjoint(
      )*Ye) + 2.7*Sqr(g1p)*(TYe*Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*Ye.adjoint()
      *Ye) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*
      Ye.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye));


   return beta_TYe;
}

} // namespace flexiblesusy
