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
 * Calculates the one-loop beta function of MuPr.
 *
 * @return one-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_MuPr_one_loop(const Susy_traces& susy_traces) const
{


   double beta_MuPr;

   beta_MuPr = -0.2*oneOver16PiSqr*MuPr*(3*Sqr(g1) + 2*Sqr(g1p) + 15*Sqr(
      g2));


   return beta_MuPr;
}

/**
 * Calculates the two-loop beta function of MuPr.
 *
 * @return two-loop beta function
 */
double NE6SSM_susy_parameters::calc_beta_MuPr_two_loop(const Susy_traces& susy_traces) const
{
   const auto QS = INPUT(QS);


   double beta_MuPr;

   beta_MuPr = 0.01*twoLoop*MuPr*(522*Power(g1,4) + 1050*Power(g2,4) +
      120*Sqr(g1p)*Sqr(g2) + 12*Sqr(g1)*(8*Sqr(g1p) + 15*Sqr(g2)) + Power(g1p,4
      )*(257 + 2*Sqr(QS)));


   return beta_MuPr;
}

} // namespace flexiblesusy
