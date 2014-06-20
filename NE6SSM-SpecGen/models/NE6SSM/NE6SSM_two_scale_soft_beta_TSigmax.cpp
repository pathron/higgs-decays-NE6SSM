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
 * Calculates the one-loop beta function of TSigmax.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_TSigmax_one_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_TSigmax;

   beta_TSigmax = oneOver16PiSqr*(6*traceAdjKappaTKappa*Sigmax + 0.2*
      MassBp*Sigmax*Sqr(g1p)*Sqr(QS) + 3*traceKappaAdjKappa*TSigmax + 9*AbsSqr(
      Sigmax)*TSigmax - 0.1*Sqr(g1p)*Sqr(QS)*TSigmax + 2*Conj(KappaPr)*(2*
      Sigmax*TKappaPr + KappaPr*TSigmax) + 2*Conj(Lambdax)*(2*Sigmax*TLambdax +
      Lambdax*TSigmax));


   return beta_TSigmax;
}

/**
 * Calculates the two-loop beta function of TSigmax.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_TSigmax_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;


   double beta_TSigmax;

   beta_TSigmax = twoLoop*(-0.04*Power(g1p,4)*MassBp*Power(QS,4)*Sigmax -
      24*traceKappaAdjKappaTKappaAdjKappa*Sigmax + 1.6*traceAdjKappaTKappa*
      Sigmax*Sqr(g1) - 1.6*MassB*traceKappaAdjKappa*Sigmax*Sqr(g1) + 3.9*
      traceAdjKappaTKappa*Sigmax*Sqr(g1p) - 3.9*MassBp*traceKappaAdjKappa*
      Sigmax*Sqr(g1p) + 32*traceAdjKappaTKappa*Sigmax*Sqr(g3) - 32*MassG*
      traceKappaAdjKappa*Sigmax*Sqr(g3) - 2.49*Power(g1p,4)*MassBp*Sigmax*Sqr(
      QS) - 0.3*traceAdjKappaTKappa*Sigmax*Sqr(g1p)*Sqr(QS) + 0.3*MassBp*
      traceKappaAdjKappa*Sigmax*Sqr(g1p)*Sqr(QS) - 32*KappaPr*Sigmax*Sqr(Conj(
      KappaPr))*TKappaPr + 0.01*Power(g1p,4)*Power(QS,4)*TSigmax - 6*
      traceKappaAdjKappaKappaAdjKappa*TSigmax + 0.8*traceKappaAdjKappa*Sqr(g1)*
      TSigmax + 1.95*traceKappaAdjKappa*Sqr(g1p)*TSigmax + 16*
      traceKappaAdjKappa*Sqr(g3)*TSigmax + 0.6225*Power(g1p,4)*Sqr(QS)*TSigmax
      - 0.15*traceKappaAdjKappa*Sqr(g1p)*Sqr(QS)*TSigmax - 8*Sqr(Conj(KappaPr))
      *Sqr(KappaPr)*TSigmax - 30*Sqr(Conj(Sigmax))*Sqr(Sigmax)*TSigmax - 4*
      Lambdax*Sqr(Conj(Lambdax))*(4*Sigmax*TLambdax + Lambdax*TSigmax) - 0.1*
      AbsSqr(Sigmax)*(120*traceAdjKappaTKappa*Sigmax + 2*MassBp*Sigmax*Sqr(g1p)
      *Sqr(QS) + 180*traceKappaAdjKappa*TSigmax - 3*Sqr(g1p)*Sqr(QS)*TSigmax +
      80*Conj(KappaPr)*(2*Sigmax*TKappaPr + 3*KappaPr*TSigmax)) - 0.1*Conj(
      Lambdax)*(2*Sigmax*(60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu +
      40*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2) + Sqr(g1p)*Sqr
      (QS))*TLambdax + Lambdax*(2*Sigmax*(60*traceAdjYdTYd + 20*traceAdjYeTYe +
      60*traceAdjYuTYu + 12*MassB*Sqr(g1) + 13*MassBp*Sqr(g1p) + 60*MassWB*Sqr
      (g2) - MassBp*Sqr(g1p)*Sqr(QS)) + (60*traceYdAdjYd + 20*traceYeAdjYe + 60
      *traceYuAdjYu + 120*AbsSqr(Sigmax) - 12*Sqr(g1) - 13*Sqr(g1p) - 60*Sqr(g2
      ) + Sqr(g1p)*Sqr(QS))*TSigmax)));


   return beta_TSigmax;
}

} // namespace flexiblesusy
