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
 * Calculates the one-loop beta function of vphi.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vphi_one_loop(const Susy_traces& susy_traces) const
{


   double beta_vphi;

   beta_vphi = -(oneOver16PiSqr*vphi*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax))
      );


   return beta_vphi;
}

/**
 * Calculates the two-loop beta function of vphi.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_vphi_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_vphi;

   beta_vphi = twoLoop*(4*vphi*AbsSqr(KappaPr)*AbsSqr(Sigmax) + 0.1*vphi*
      AbsSqr(Sigmax)*(30*traceKappaAdjKappa + 20*AbsSqr(Lambdax) + 20*AbsSqr(
      Sigmax) - Sqr(g1p)*Sqr(QS)) + 8*vphi*Sqr(Conj(KappaPr))*Sqr(KappaPr));


   return beta_vphi;
}

} // namespace flexiblesusy
