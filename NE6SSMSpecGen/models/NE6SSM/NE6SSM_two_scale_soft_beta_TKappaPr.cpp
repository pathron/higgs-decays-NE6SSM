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

// File generated at Fri 13 Jun 2014 15:36:37

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TKappaPr.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_TKappaPr_one_loop(const Soft_traces& soft_traces) const
{


   double beta_TKappaPr;

   beta_TKappaPr = 3*oneOver16PiSqr*(6*AbsSqr(KappaPr)*TKappaPr + Conj(
      Sigmax)*(Sigmax*TKappaPr + 2*KappaPr*TSigmax));


   return beta_TKappaPr;
}

/**
 * Calculates the two-loop beta function of TKappaPr.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_TKappaPr_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_TKappaPr;

   beta_TKappaPr = -0.3*twoLoop*(400*Sqr(Conj(KappaPr))*Sqr(KappaPr)*
      TKappaPr + 20*Sigmax*Sqr(Conj(Sigmax))*(Sigmax*TKappaPr + 4*KappaPr*
      TSigmax) + Conj(Sigmax)*(Sigmax*(30*traceKappaAdjKappa + 120*AbsSqr(
      KappaPr) + 20*AbsSqr(Lambdax) - Sqr(g1p)*Sqr(QS))*TKappaPr + 2*KappaPr*(
      Sigmax*(30*traceAdjKappaTKappa + MassBp*Sqr(g1p)*Sqr(QS)) + (30*
      traceKappaAdjKappa + 40*AbsSqr(KappaPr) - Sqr(g1p)*Sqr(QS))*TSigmax + 20*
      Conj(Lambdax)*(Sigmax*TLambdax + Lambdax*TSigmax))));


   return beta_TKappaPr;
}

} // namespace flexiblesusy
