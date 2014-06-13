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

// File generated at Fri 13 Jun 2014 15:36:47

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassWB.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_MassWB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = 4*MassWB*oneOver16PiSqr*Sqr(g2);


   return beta_MassWB;
}

/**
 * Calculates the two-loop beta function of MassWB.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_MassWB_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_MassWB;

   beta_MassWB = 0.2*twoLoop*Sqr(g2)*(60*traceAdjYdTYd + 20*traceAdjYeTYe
      + 60*traceAdjYuTYu - 60*MassWB*traceYdAdjYd - 20*MassWB*traceYeAdjYe -
      60*MassWB*traceYuAdjYu + 24*MassB*Sqr(g1) + 24*MassWB*Sqr(g1) + 21*MassBp
      *Sqr(g1p) + 21*MassWB*Sqr(g1p) + 640*MassWB*Sqr(g2) + 240*MassG*Sqr(g3) +
      240*MassWB*Sqr(g3) - 20*Conj(Lambdax)*(MassWB*Lambdax - TLambdax));


   return beta_MassWB;
}

} // namespace flexiblesusy
