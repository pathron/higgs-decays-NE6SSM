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
 * Calculates the one-loop beta function of KappaPr.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_KappaPr_one_loop(const Susy_traces& susy_traces) const
{


   double beta_KappaPr;

   beta_KappaPr = 3*oneOver16PiSqr*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax))*
      KappaPr;


   return beta_KappaPr;
}

/**
 * Calculates the two-loop beta function of KappaPr.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_KappaPr_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_KappaPr;

   beta_KappaPr = -0.3*twoLoop*KappaPr*(40*AbsSqr(KappaPr)*AbsSqr(Sigmax)
      + AbsSqr(Sigmax)*(30*traceKappaAdjKappa + 20*AbsSqr(Lambdax) + 20*AbsSqr
      (Sigmax) - Sqr(g1p)*Sqr(QS)) + 80*Sqr(Conj(KappaPr))*Sqr(KappaPr));


   return beta_KappaPr;
}

} // namespace flexiblesusy
