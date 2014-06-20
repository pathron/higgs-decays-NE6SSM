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

// File generated at Fri 20 Jun 2014 21:37:54

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassB.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_MassB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = 16.8*MassB*oneOver16PiSqr*Sqr(g1);


   return beta_MassB;
}

/**
 * Calculates the two-loop beta function of MassB.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_MassB_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_MassB;

   beta_MassB = 0.04*twoLoop*Sqr(g1)*(40*traceAdjKappaTKappa + 140*
      traceAdjYdTYd + 180*traceAdjYeTYe + 260*traceAdjYuTYu - 40*MassB*
      traceKappaAdjKappa - 140*MassB*traceYdAdjYd - 180*MassB*traceYeAdjYe -
      260*MassB*traceYuAdjYu + 864*MassB*Sqr(g1) + 123*MassB*Sqr(g1p) + 123*
      MassBp*Sqr(g1p) + 360*MassB*Sqr(g2) + 360*MassWB*Sqr(g2) + 1200*MassB*Sqr
      (g3) + 1200*MassG*Sqr(g3) - 60*Conj(Lambdax)*(MassB*Lambdax - TLambdax));


   return beta_MassB;
}

} // namespace flexiblesusy
