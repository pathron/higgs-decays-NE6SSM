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
 * Calculates the one-loop beta function of MassBp.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_MassBp_one_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);


   double beta_MassBp;

   beta_MassBp = 0.05*MassBp*oneOver16PiSqr*Sqr(g1p)*(249 + 2*Sqr(QS));


   return beta_MassBp;
}

/**
 * Calculates the two-loop beta function of MassBp.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_MassBp_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_MassBp;

   beta_MassBp = 0.01*twoLoop*Sqr(g1p)*(390*traceAdjKappaTKappa + 840*
      traceAdjYdTYd + 280*traceAdjYeTYe + 360*traceAdjYuTYu - 390*MassBp*
      traceKappaAdjKappa - 840*MassBp*traceYdAdjYd - 280*MassBp*traceYeAdjYe -
      360*MassBp*traceYuAdjYu + 492*MassB*Sqr(g1) + 492*MassBp*Sqr(g1) + 1401*
      MassBp*Sqr(g1p) + 2*MassBp*Power(QS,4)*Sqr(g1p) + 1260*MassBp*Sqr(g2) +
      1260*MassWB*Sqr(g2) + 4800*MassBp*Sqr(g3) + 4800*MassG*Sqr(g3) + 30*
      traceAdjKappaTKappa*Sqr(QS) - 30*MassBp*traceKappaAdjKappa*Sqr(QS) - 20*
      Conj(Lambdax)*(13 + Sqr(QS))*(MassBp*Lambdax - TLambdax) - 20*Conj(Sigmax
      )*Sqr(QS)*(MassBp*Sigmax - TSigmax));


   return beta_MassBp;
}

} // namespace flexiblesusy
