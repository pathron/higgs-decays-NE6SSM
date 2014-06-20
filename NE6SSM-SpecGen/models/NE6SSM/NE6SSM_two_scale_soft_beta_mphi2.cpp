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
 * Calculates the one-loop beta function of mphi2.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_mphi2_one_loop(const Soft_traces& soft_traces) const
{


   double beta_mphi2;

   beta_mphi2 = 2*oneOver16PiSqr*(6*mphi2*AbsSqr(KappaPr) + (mphi2 + ms2
      + msbar2)*AbsSqr(Sigmax) + 2*AbsSqr(TKappaPr) + AbsSqr(TSigmax));


   return beta_mphi2;
}

/**
 * Calculates the two-loop beta function of mphi2.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_mphi2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;


   double beta_mphi2;

   beta_mphi2 = twoLoop*(-96*mphi2*Sqr(Conj(KappaPr))*Sqr(KappaPr) - 8*(
      mphi2 + ms2 + msbar2)*Sqr(Conj(Sigmax))*Sqr(Sigmax) + 0.2*Conj(Sigmax)*(
      -30*traceconjTKappaTpTKappa*Sigmax - 30*mphi2*traceKappaAdjKappa*Sigmax -
      60*ms2*traceKappaAdjKappa*Sigmax - 30*msbar2*traceKappaAdjKappa*Sigmax -
      30*traceKappaAdjKappaconjmDx2*Sigmax - 30*traceKappaconjmDxbar2AdjKappa*
      Sigmax - 20*(mHd2 + mHu2 + mphi2 + 2*ms2 + msbar2)*AbsSqr(Lambdax)*Sigmax
      - 40*AbsSqr(TKappaPr)*Sigmax - 20*AbsSqr(TLambdax)*Sigmax - 80*AbsSqr(
      TSigmax)*Sigmax + mphi2*Sigmax*Sqr(g1p)*Sqr(QS) + ms2*Sigmax*Sqr(g1p)*Sqr
      (QS) + msbar2*Sigmax*Sqr(g1p)*Sqr(QS) + Conj(MassBp)*Sqr(g1p)*Sqr(QS)*(2*
      MassBp*Sigmax - TSigmax) - 30*traceconjTKappaTpKappa*TSigmax - 40*Conj(
      TKappaPr)*KappaPr*TSigmax - 20*Conj(TLambdax)*Lambdax*TSigmax) - 8*Conj(
      KappaPr)*((4*mphi2 + ms2 + msbar2)*AbsSqr(Sigmax)*KappaPr + 8*AbsSqr(
      TKappaPr)*KappaPr + Conj(TSigmax)*(Sigmax*TKappaPr + KappaPr*TSigmax)) -
      0.2*Conj(TSigmax)*(30*traceAdjKappaTKappa*Sigmax + MassBp*Sigmax*Sqr(g1p)
      *Sqr(QS) + 30*traceKappaAdjKappa*TSigmax - Sqr(g1p)*Sqr(QS)*TSigmax + 20*
      Conj(Lambdax)*(Sigmax*TLambdax + Lambdax*TSigmax)));


   return beta_mphi2;
}

} // namespace flexiblesusy
