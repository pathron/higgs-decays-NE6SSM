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

// File generated at Fri 20 Jun 2014 21:37:46

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of BMuPhi.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_BMuPhi_one_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPhi;

   beta_BMuPhi = oneOver16PiSqr*(8*AbsSqr(KappaPr)*BMuPhi + 2*AbsSqr(
      Sigmax)*BMuPhi + 8*MuPhi*Conj(KappaPr)*TKappaPr + 4*MuPhi*Conj(Sigmax)*
      TSigmax);


   return beta_BMuPhi;
}

/**
 * Calculates the two-loop beta function of BMuPhi.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_BMuPhi_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_BMuPhi;

   beta_BMuPhi = twoLoop*(-0.2*BMuPhi*(80*AbsSqr(KappaPr)*AbsSqr(Sigmax)
      + AbsSqr(Sigmax)*(30*traceKappaAdjKappa + 20*AbsSqr(Lambdax) + 20*AbsSqr(
      Sigmax) - Sqr(g1p)*Sqr(QS)) + 160*Sqr(Conj(KappaPr))*Sqr(KappaPr)) - 0.4*
      MuPhi*(200*KappaPr*Sqr(Conj(KappaPr))*TKappaPr + 40*Sigmax*Sqr(Conj(
      Sigmax))*TSigmax + Conj(Sigmax)*(30*traceAdjKappaTKappa*Sigmax + MassBp*
      Sigmax*Sqr(g1p)*Sqr(QS) + 30*traceKappaAdjKappa*TSigmax - Sqr(g1p)*Sqr(QS
      )*TSigmax + 20*Conj(KappaPr)*(2*Sigmax*TKappaPr + 3*KappaPr*TSigmax) + 20
      *Conj(Lambdax)*(Sigmax*TLambdax + Lambdax*TSigmax))));


   return beta_BMuPhi;
}

} // namespace flexiblesusy
