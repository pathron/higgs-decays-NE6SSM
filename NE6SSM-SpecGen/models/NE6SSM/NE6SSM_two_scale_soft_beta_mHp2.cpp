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

// File generated at Fri 20 Jun 2014 21:37:53

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHp2.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_mHp2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHp2;

   beta_mHp2 = 0.2*oneOver16PiSqr*(-3.872983346207417*g1*Tr11 +
      3.1622776601683795*g1p*Tr14 - 6*AbsSqr(MassB)*Sqr(g1) - 4*AbsSqr(MassBp)*
      Sqr(g1p) - 30*AbsSqr(MassWB)*Sqr(g2));


   return beta_mHp2;
}

/**
 * Calculates the two-loop beta function of mHp2.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_mHp2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHp2;

   beta_mHp2 = 0.02*twoLoop*(6*Conj(MassB)*Sqr(g1)*(261*MassB*Sqr(g1) + 8
      *(2*MassB + MassBp)*Sqr(g1p) + 15*(2*MassB + MassWB)*Sqr(g2)) + 10*(30*
      Power(g2,4)*Tr22 - 4.898979485566356*g1*g1p*Tr2U114 - 4.898979485566356*
      g1*g1p*Tr2U141 - 15.491933384829668*g1*Tr31 + 12.649110640673518*g1p*Tr34
      + 6*Tr2U111*Sqr(g1) + 4*Tr2U144*Sqr(g1p) + 3*Conj(MassWB)*Sqr(g2)*(3*(
      MassB + 2*MassWB)*Sqr(g1) + 2*(MassBp + 2*MassWB)*Sqr(g1p) + 85*MassWB*
      Sqr(g2))) + 3*Conj(MassBp)*Sqr(g1p)*(16*(MassB + 2*MassBp)*Sqr(g1) + 20*(
      2*MassBp + MassWB)*Sqr(g2) + MassBp*Sqr(g1p)*(257 + 2*Sqr(QS))));


   return beta_mHp2;
}

} // namespace flexiblesusy
