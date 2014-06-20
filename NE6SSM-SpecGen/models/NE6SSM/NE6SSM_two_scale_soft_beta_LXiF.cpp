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
 * Calculates the one-loop beta function of LXiF.
 *
 * @return one-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_LXiF_one_loop(const Soft_traces& soft_traces) const
{


   double beta_LXiF;

   beta_LXiF = oneOver16PiSqr*(2*MuPhi*BMuPhi*Conj(KappaPr) + 4*mphi2*
      Conj(MuPhi)*KappaPr + 2*AbsSqr(KappaPr)*LXiF + AbsSqr(Sigmax)*LXiF + 2*
      Conj(BMuPhi)*TKappaPr + 4*Conj(KappaPr)*XiF*TKappaPr + 2*Conj(Sigmax)*XiF
      *TSigmax);


   return beta_LXiF;
}

/**
 * Calculates the two-loop beta function of LXiF.
 *
 * @return two-loop beta function
 */
double NE6SSM_soft_parameters::calc_beta_LXiF_two_loop(const Soft_traces& soft_traces) const
{
   const auto QS = INPUT(QS);
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;


   double beta_LXiF;

   beta_LXiF = twoLoop*(-4*MuPhi*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax))*
      BMuPhi*Conj(KappaPr) - 6*traceAdjKappaTKappa*AbsSqr(Sigmax)*XiF - 3*
      traceKappaAdjKappa*AbsSqr(Sigmax)*LXiF - 4*AbsSqr(KappaPr)*AbsSqr(Sigmax)
      *LXiF - 2*AbsSqr(Lambdax)*AbsSqr(Sigmax)*LXiF - 0.2*MassBp*AbsSqr(Sigmax)
      *XiF*Sqr(g1p)*Sqr(QS) + 0.1*AbsSqr(Sigmax)*LXiF*Sqr(g1p)*Sqr(QS) - 8*LXiF
      *Sqr(Conj(KappaPr))*Sqr(KappaPr) - 2*LXiF*Sqr(Conj(Sigmax))*Sqr(Sigmax) -
      16*AbsSqr(KappaPr)*Conj(BMuPhi)*TKappaPr - 4*AbsSqr(Sigmax)*Conj(BMuPhi)
      *TKappaPr - 8*AbsSqr(Sigmax)*Conj(KappaPr)*XiF*TKappaPr - 32*KappaPr*XiF*
      Sqr(Conj(KappaPr))*TKappaPr - 8*Sqr(MuPhi)*Sqr(Conj(KappaPr))*TKappaPr -
      4*Conj(MuPhi)*((3*mphi2 + ms2 + msbar2)*AbsSqr(Sigmax)*KappaPr + 4*AbsSqr
      (TKappaPr)*KappaPr + AbsSqr(TSigmax)*KappaPr + 10*mphi2*Conj(KappaPr)*Sqr
      (KappaPr) + Conj(TSigmax)*Sigmax*TKappaPr) - 4*AbsSqr(Sigmax)*Conj(
      Lambdax)*XiF*TLambdax - 4*Conj(BMuPhi)*Conj(Sigmax)*KappaPr*TSigmax - 6*
      traceKappaAdjKappa*Conj(Sigmax)*XiF*TSigmax - 8*AbsSqr(KappaPr)*Conj(
      Sigmax)*XiF*TSigmax - 4*AbsSqr(Lambdax)*Conj(Sigmax)*XiF*TSigmax - 4*Conj
      (KappaPr)*Conj(Sigmax)*Sqr(MuPhi)*TSigmax + 0.2*Conj(Sigmax)*XiF*Sqr(g1p)
      *Sqr(QS)*TSigmax - 8*XiF*Sigmax*Sqr(Conj(Sigmax))*TSigmax);


   return beta_LXiF;
}

} // namespace flexiblesusy
