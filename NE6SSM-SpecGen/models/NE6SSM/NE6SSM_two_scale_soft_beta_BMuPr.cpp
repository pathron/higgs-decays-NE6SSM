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
 * Calculates the one-loop beta function of BMuPr.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_BMuPr_one_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = 0.2*oneOver16PiSqr*(-(BMuPr*(3*Sqr(g1) + 2*Sqr(g1p) + 15*
      Sqr(g2))) + 2*MuPr*(3*MassB*Sqr(g1) + 2*MassBp*Sqr(g1p) + 15*MassWB*Sqr(
      g2)));


   return beta_BMuPr;
}

/**
 * Calculates the two-loop beta function of BMuPr.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_BMuPr_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);


   double beta_BMuPr;

   beta_BMuPr = 0.01*twoLoop*(BMuPr*(522*Power(g1,4) + 1050*Power(g2,4) +
      120*Sqr(g1p)*Sqr(g2) + 12*Sqr(g1)*(2*Sqr(g1p) + 15*Sqr(g2)) + Power(g1p,
      4)*(257 + 2*Sqr(QS))) - 4*MuPr*(522*Power(g1,4)*MassB + 1050*Power(g2,4)*
      MassWB + 60*(MassBp + MassWB)*Sqr(g1p)*Sqr(g2) + 6*Sqr(g1)*(2*(MassB +
      MassBp)*Sqr(g1p) + 15*(MassB + MassWB)*Sqr(g2)) + Power(g1p,4)*MassBp*(
      257 + 2*Sqr(QS))));


   return beta_BMuPr;
}

} // namespace flexiblesusy
