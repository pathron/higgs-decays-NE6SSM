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

// File generated at Fri 13 Jun 2014 15:35:55

#include "NE6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME NE6SSM_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

NE6SSM_susy_parameters::NE6SSM_susy_parameters(const NE6SSM_input_parameters& input_)
   : Beta_function()
   , Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero(
   )), KappaPr(0), Sigmax(0), Kappa(Eigen::Matrix<double,3,3>::Zero()), Lambdax
   (0), Yu(Eigen::Matrix<double,3,3>::Zero()), MuPr(0), MuPhi(0), XiF(0), g1(0)
   , g2(0), g3(0), g1p(0), vd(0), vu(0), vs(0), vsb(0), vphi(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

NE6SSM_susy_parameters::NE6SSM_susy_parameters(
   double scale_, double loops_, double thresholds_,
   const NE6SSM_input_parameters& input_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , double KappaPr_, double Sigmax_, const Eigen::Matrix<double,3,3>& Kappa_,
   double Lambdax_, const Eigen::Matrix<double,3,3>& Yu_, double MuPr_, double
   MuPhi_, double XiF_, double g1_, double g2_, double g3_, double g1p_, double
   vd_, double vu_, double vs_, double vsb_, double vphi_

)
   : Beta_function()
   , Yd(Yd_), Ye(Ye_), KappaPr(KappaPr_), Sigmax(Sigmax_), Kappa(Kappa_),
   Lambdax(Lambdax_), Yu(Yu_), MuPr(MuPr_), MuPhi(MuPhi_), XiF(XiF_), g1(g1_),
   g2(g2_), g3(g3_), g1p(g1p_), vd(vd_), vu(vu_), vs(vs_), vsb(vsb_), vphi(
   vphi_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd NE6SSM_susy_parameters::beta() const
{
   return calc_beta().get();
}

NE6SSM_susy_parameters NE6SSM_susy_parameters::calc_beta() const
{
   Susy_traces susy_traces;
   calc_susy_traces(susy_traces);

   Eigen::Matrix<double,3,3> beta_Yd(calc_beta_Yd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Ye(calc_beta_Ye_one_loop(TRACE_STRUCT));
   double beta_KappaPr(calc_beta_KappaPr_one_loop(TRACE_STRUCT));
   double beta_Sigmax(calc_beta_Sigmax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Kappa(calc_beta_Kappa_one_loop(TRACE_STRUCT))
      ;
   double beta_Lambdax(calc_beta_Lambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Yu(calc_beta_Yu_one_loop(TRACE_STRUCT));
   double beta_MuPr(calc_beta_MuPr_one_loop(TRACE_STRUCT));
   double beta_MuPhi(calc_beta_MuPhi_one_loop(TRACE_STRUCT));
   double beta_XiF(calc_beta_XiF_one_loop(TRACE_STRUCT));
   double beta_g1(calc_beta_g1_one_loop(TRACE_STRUCT));
   double beta_g2(calc_beta_g2_one_loop(TRACE_STRUCT));
   double beta_g3(calc_beta_g3_one_loop(TRACE_STRUCT));
   double beta_g1p(calc_beta_g1p_one_loop(TRACE_STRUCT));
   double beta_vd(calc_beta_vd_one_loop(TRACE_STRUCT));
   double beta_vu(calc_beta_vu_one_loop(TRACE_STRUCT));
   double beta_vs(calc_beta_vs_one_loop(TRACE_STRUCT));
   double beta_vsb(calc_beta_vsb_one_loop(TRACE_STRUCT));
   double beta_vphi(calc_beta_vphi_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
      beta_KappaPr += calc_beta_KappaPr_two_loop(TRACE_STRUCT);
      beta_Sigmax += calc_beta_Sigmax_two_loop(TRACE_STRUCT);
      beta_Kappa += calc_beta_Kappa_two_loop(TRACE_STRUCT);
      beta_Lambdax += calc_beta_Lambdax_two_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
      beta_MuPr += calc_beta_MuPr_two_loop(TRACE_STRUCT);
      beta_MuPhi += calc_beta_MuPhi_two_loop(TRACE_STRUCT);
      beta_XiF += calc_beta_XiF_two_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
      beta_g1p += calc_beta_g1p_two_loop(TRACE_STRUCT);
      beta_vd += calc_beta_vd_two_loop(TRACE_STRUCT);
      beta_vu += calc_beta_vu_two_loop(TRACE_STRUCT);
      beta_vs += calc_beta_vs_two_loop(TRACE_STRUCT);
      beta_vsb += calc_beta_vsb_two_loop(TRACE_STRUCT);
      beta_vphi += calc_beta_vphi_two_loop(TRACE_STRUCT);

   }


   return NE6SSM_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_Yd, beta_Ye, beta_KappaPr, beta_Sigmax, beta_Kappa, beta_Lambdax, beta_Yu, beta_MuPr, beta_MuPhi, beta_XiF, beta_g1, beta_g2, beta_g3, beta_g1p, beta_vd, beta_vu, beta_vs, beta_vsb, beta_vphi);
}

void NE6SSM_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   KappaPr = 0.;
   Sigmax = 0.;
   Kappa = Eigen::Matrix<double,3,3>::Zero();
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   MuPr = 0.;
   MuPhi = 0.;
   XiF = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   g1p = 0.;
   vd = 0.;
   vu = 0.;
   vs = 0.;
   vsb = 0.;
   vphi = 0.;

}

Eigen::Matrix<double,3,3> CLASSNAME::get_SqSq() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(Yd.adjoint()*Yd + Yu.adjoint()*Yu -
      0.016666666666666666*(2*Sqr(g1) + 3*Sqr(g1p) + 90*Sqr(g2) + 160*Sqr(g3))*
      UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(0.2*(-10*(Yd.adjoint()*Yd*Yd.adjoint()*Yd +
         Yu.adjoint()*Yu*Yu.adjoint()*Yu) + Yd.adjoint()*Yd*(-5*AbsSqr(Lambdax)
         + 2*Sqr(g1) + 3*Sqr(g1p) - 15*(Yd*Yd.adjoint()).trace() - 5*(Ye*
         Ye.adjoint()).trace()) + Yu.adjoint()*Yu*(-5*AbsSqr(Lambdax) + 4*Sqr(
         g1) + Sqr(g1p) - 15*(Yu*Yu.adjoint()).trace())) +
         0.0001388888888888889*(2024*Power(g1,4) - 8*Sqr(g1)*(51*Sqr(g1p) - 90*
         Sqr(g2) - 160*Sqr(g3)) + 120*Sqr(g1p)*(9*Sqr(g2) + 16*Sqr(g3)) + 200*(
         189*Power(g2,4) + 256*Power(g3,4) + 288*Sqr(g2)*Sqr(g3)) + 9*Power(g1p
         ,4)*(251 + 2*Sqr(QS)))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SlSl() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(Ye.adjoint()*Ye - 0.1*(3*Sqr(g1) + 2*Sqr(g1p
      ) + 15*Sqr(g2))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += 0.005*twoLoop*(20*(-20*(Ye.adjoint()*Ye*Ye.adjoint()*
         Ye) + Ye.adjoint()*Ye*(-10*AbsSqr(Lambdax) + 12*Sqr(g1) + 3*Sqr(g1p) -
         30*(Yd*Yd.adjoint()).trace() - 10*(Ye*Ye.adjoint()).trace())) + (522*
         Power(g1,4) + 1050*Power(g2,4) + 120*Sqr(g1p)*Sqr(g2) + 12*Sqr(g1)*(8*
         Sqr(g1p) + 15*Sqr(g2)) + Power(g1p,4)*(257 + 2*Sqr(QS)))*UNITMATRIX(3)
         );
   }

   return anomDim;
}

double CLASSNAME::get_SHdSHd() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 0.45*Sqr(g1p
      ) - 1.5*Sqr(g2) + 3*(Yd*Yd.adjoint()).trace() + (Ye*Ye.adjoint()).trace()
      );

   if (get_loops() > 1) {
      anomDim += twoLoop*(-3*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.05*
         AbsSqr(Lambdax)*(-20*AbsSqr(Sigmax) - 5*Sqr(g1p) + Sqr(g1p)*Sqr(QS) -
         60*(Yu*Yu.adjoint()).trace() - 60*(Kappa*(Kappa).adjoint()).trace()) +
         0.00125*(-160*(2*Sqr(g1) + 3*Sqr(g1p) - 80*Sqr(g3))*(Yd*Yd.adjoint())
         .trace() + 160*(6*Sqr(g1) - Sqr(g1p))*(Ye*Ye.adjoint()).trace() + 3*(
         696*Power(g1,4) + 801*Power(g1p,4) + 1400*Power(g2,4) - 72*Sqr(g1)*Sqr
         (g1p) + 240*Sqr(g1)*Sqr(g2) + 360*Sqr(g1p)*Sqr(g2) + 6*Power(g1p,4)*
         Sqr(QS) - 2400*(Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace() - 800*(Yd*
         Yu.adjoint()*Yu*Yd.adjoint()).trace() - 800*(Ye*Ye.adjoint()*Ye*
         Ye.adjoint()).trace())));
   }

   return anomDim;
}

double CLASSNAME::get_SHuSHu() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(AbsSqr(Lambdax) - 0.3*Sqr(g1) - 0.2*Sqr(g1p)
      - 1.5*Sqr(g2) + 3*(Yu*Yu.adjoint()).trace());

   if (get_loops() > 1) {
      anomDim += 0.005*twoLoop*(522*Power(g1,4) + 257*Power(g1p,4) +
         1050*Power(g2,4) + 96*Sqr(g1)*Sqr(g1p) + 180*Sqr(g1)*Sqr(g2) + 120*Sqr
         (g1p)*Sqr(g2) + 2*Power(g1p,4)*Sqr(QS) - 600*Sqr(Conj(Lambdax))*Sqr(
         Lambdax) + 20*(8*Sqr(g1) - 3*Sqr(g1p) + 160*Sqr(g3))*(Yu*Yu.adjoint())
         .trace() + 10*AbsSqr(Lambdax)*(-20*AbsSqr(Sigmax) + 5*Sqr(g1p) + Sqr(
         g1p)*Sqr(QS) - 60*(Yd*Yd.adjoint()).trace() - 20*(Ye*Ye.adjoint())
         .trace() - 60*(Kappa*(Kappa).adjoint()).trace()) - 600*(Yd*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace() - 1800*(Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SdRSdR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = 0.06666666666666667*oneOver16PiSqr*(30*(Yd.conjugate()*
      Yd.transpose()) - (2*Sqr(g1) + 3*Sqr(g1p) + 40*Sqr(g3))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-2*(Yd.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yd.transpose() + Yd.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yd.transpose()) + Yd.conjugate()*Yd.transpose()*(-2*
         AbsSqr(Lambdax) + 0.4*Sqr(g1) + 0.6*Sqr(g1p) + 6*Sqr(g2) - 6*(Yd*
         Yd.adjoint()).trace() - 2*(Ye*Ye.adjoint()).trace()) +
         0.0005555555555555556*(2048*Power(g1,4) + 12800*Power(g3,4) - 16*Sqr(
         g1)*(21*Sqr(g1p) - 80*Sqr(g3)) + 1920*Sqr(g1p)*Sqr(g3) + 9*Power(g1p,4
         )*(257 + 2*Sqr(QS)))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SuRSuR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(2*(Yu.conjugate()*Yu.transpose()) -
      0.016666666666666666*(32*Sqr(g1) + 3*Sqr(g1p) + 160*Sqr(g3))*UNITMATRIX(3
      ));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-0.4*(5*(Yu.conjugate()*Yd.transpose()*
         Yd.conjugate()*Yu.transpose() + Yu.conjugate()*Yu.transpose()*
         Yu.conjugate()*Yu.transpose()) + Yu.conjugate()*Yu.transpose()*(5*
         AbsSqr(Lambdax) + Sqr(g1) - Sqr(g1p) - 15*Sqr(g2) + 15*(Yu*Yu.adjoint(
         )).trace())) + 0.0001388888888888889*(34304*Power(g1,4) + 51200*Power(
         g3,4) + 1920*Sqr(g1p)*Sqr(g3) + 64*Sqr(g1)*(33*Sqr(g1p) + 320*Sqr(g3))
         + 9*Power(g1p,4)*(251 + 2*Sqr(QS)))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SeRSeR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(2*(Ye.conjugate()*Ye.transpose()) - 0.05*(24
      *Sqr(g1) + Sqr(g1p))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-0.4*(5*(Ye.conjugate()*Ye.transpose()*
         Ye.conjugate()*Ye.transpose()) + Ye.conjugate()*Ye.transpose()*(5*
         AbsSqr(Lambdax) + 3*Sqr(g1) - 3*Sqr(g1p) - 15*Sqr(g2) + 15*(Yd*
         Yd.adjoint()).trace() + 5*(Ye*Ye.adjoint()).trace())) + 0.00125*(9216*
         Power(g1,4) - 192*Sqr(g1)*Sqr(g1p) + Power(g1p,4)*(251 + 2*Sqr(QS)))*
         UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SsRSsR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(2*AbsSqr(Lambdax) + AbsSqr(Sigmax) - 0.05*
      Sqr(g1p)*Sqr(QS) + 3*(Kappa*(Kappa).adjoint()).trace());

   if (get_loops() > 1) {
      anomDim += twoLoop*(0.005*Power(g1p,4)*Power(QS,4) - 2*AbsSqr(
         KappaPr)*AbsSqr(Sigmax) + 0.31125*Power(g1p,4)*Sqr(QS) - 4*Sqr(Conj(
         Lambdax))*Sqr(Lambdax) - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax) + 0.1*AbsSqr(
         Lambdax)*(12*Sqr(g1) + 13*Sqr(g1p) + 60*Sqr(g2) - Sqr(g1p)*Sqr(QS) -
         60*(Yd*Yd.adjoint()).trace() - 20*(Ye*Ye.adjoint()).trace() - 60*(Yu*
         Yu.adjoint()).trace()) + 0.8*Sqr(g1)*(Kappa*(Kappa).adjoint()).trace()
         + 1.95*Sqr(g1p)*(Kappa*(Kappa).adjoint()).trace() + 16*Sqr(g3)*(Kappa
         *(Kappa).adjoint()).trace() - 0.15*Sqr(g1p)*Sqr(QS)*(Kappa*(Kappa)
         .adjoint()).trace() - 6*(Kappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint
         ()).trace());
   }

   return anomDim;
}

double CLASSNAME::get_SsbarRSsbarR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(AbsSqr(Sigmax) - 0.05*Sqr(g1p)*Sqr(QS));

   if (get_loops() > 1) {
      anomDim += twoLoop*(0.005*Power(g1p,4)*Power(QS,4) - 2*AbsSqr(
         KappaPr)*AbsSqr(Sigmax) - 2*AbsSqr(Lambdax)*AbsSqr(Sigmax) + 0.31125*
         Power(g1p,4)*Sqr(QS) - 2*Sqr(Conj(Sigmax))*Sqr(Sigmax) - 3*AbsSqr(
         Sigmax)*(Kappa*(Kappa).adjoint()).trace());
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxLSDxL() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(Kappa.conjugate()*(Kappa).transpose() -
      0.06666666666666667*(2*Sqr(g1) + 3*Sqr(g1p) + 40*Sqr(g3))*UNITMATRIX(3));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-(Kappa.conjugate()*(Kappa).transpose()*
         Kappa.conjugate()*(Kappa).transpose()) + Kappa.conjugate()*(Kappa)
         .transpose()*(-2*AbsSqr(Lambdax) - AbsSqr(Sigmax) + 0.25*Sqr(g1p) +
         0.05*Sqr(g1p)*Sqr(QS) - 3*(Kappa*(Kappa).adjoint()).trace()) +
         0.0005555555555555556*(2048*Power(g1,4) + 12800*Power(g3,4) - 16*Sqr(
         g1)*(21*Sqr(g1p) - 80*Sqr(g3)) + 1920*Sqr(g1p)*Sqr(g3) + 9*Power(g1p,4
         )*(257 + 2*Sqr(QS)))*UNITMATRIX(3));
   }

   return anomDim;
}

Eigen::Matrix<double,3,3> CLASSNAME::get_SDxbarRSDxbarR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   Eigen::Matrix<double,3,3> anomDim;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*((Kappa).adjoint()*Kappa -
      0.016666666666666666*(8*Sqr(g1) + 27*Sqr(g1p) + 160*Sqr(g3))*UNITMATRIX(3
      ));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-((Kappa).adjoint()*Kappa*(Kappa).adjoint()*
         Kappa) + (Kappa).adjoint()*Kappa*(-2*AbsSqr(Lambdax) - AbsSqr(Sigmax)
         - 0.25*Sqr(g1p) + 0.05*Sqr(g1p)*Sqr(QS) - 3*(Kappa*(Kappa).adjoint())
         .trace()) + 0.0001388888888888889*(8192*Power(g1,4) + 51200*Power(g3,4
         ) + 17280*Sqr(g1p)*Sqr(g3) + 128*Sqr(g1)*(27*Sqr(g1p) + 40*Sqr(g3)) +
         81*Power(g1p,4)*(267 + 2*Sqr(QS)))*UNITMATRIX(3));
   }

   return anomDim;
}

double CLASSNAME::get_SHpSHp() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = 0.1*oneOver16PiSqr*(-3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2));

   if (get_loops() > 1) {
      anomDim += 0.005*twoLoop*(522*Power(g1,4) + 1050*Power(g2,4) +
         120*Sqr(g1p)*Sqr(g2) + 12*Sqr(g1)*(8*Sqr(g1p) + 15*Sqr(g2)) + Power(
         g1p,4)*(257 + 2*Sqr(QS)));
   }

   return anomDim;
}

double CLASSNAME::get_SHpbarSHpbar() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = 0.1*oneOver16PiSqr*(-3*Sqr(g1) - 2*Sqr(g1p) - 15*Sqr(g2));

   if (get_loops() > 1) {
      anomDim += 0.005*twoLoop*(522*Power(g1,4) + 1050*Power(g2,4) +
         120*Sqr(g1p)*Sqr(g2) + 12*Sqr(g1)*(8*Sqr(g1p) + 15*Sqr(g2)) + Power(
         g1p,4)*(257 + 2*Sqr(QS)));
   }

   return anomDim;
}

double CLASSNAME::get_SphiRSphiR() const
{
   const double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
   double anomDim = 0;
   const auto QS = INPUT(QS);

   anomDim = oneOver16PiSqr*(2*AbsSqr(KappaPr) + AbsSqr(Sigmax));

   if (get_loops() > 1) {
      anomDim += twoLoop*(-4*AbsSqr(KappaPr)*AbsSqr(Sigmax) - 8*Sqr(
         Conj(KappaPr))*Sqr(KappaPr) - 0.1*AbsSqr(Sigmax)*(20*AbsSqr(Lambdax) +
         20*AbsSqr(Sigmax) - Sqr(g1p)*Sqr(QS) + 30*(Kappa*(Kappa).adjoint())
         .trace()));
   }

   return anomDim;
}


const Eigen::ArrayXd NE6SSM_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Yd(0,0);
   pars(1) = Yd(0,1);
   pars(2) = Yd(0,2);
   pars(3) = Yd(1,0);
   pars(4) = Yd(1,1);
   pars(5) = Yd(1,2);
   pars(6) = Yd(2,0);
   pars(7) = Yd(2,1);
   pars(8) = Yd(2,2);
   pars(9) = Ye(0,0);
   pars(10) = Ye(0,1);
   pars(11) = Ye(0,2);
   pars(12) = Ye(1,0);
   pars(13) = Ye(1,1);
   pars(14) = Ye(1,2);
   pars(15) = Ye(2,0);
   pars(16) = Ye(2,1);
   pars(17) = Ye(2,2);
   pars(18) = KappaPr;
   pars(19) = Sigmax;
   pars(20) = Kappa(0,0);
   pars(21) = Kappa(0,1);
   pars(22) = Kappa(0,2);
   pars(23) = Kappa(1,0);
   pars(24) = Kappa(1,1);
   pars(25) = Kappa(1,2);
   pars(26) = Kappa(2,0);
   pars(27) = Kappa(2,1);
   pars(28) = Kappa(2,2);
   pars(29) = Lambdax;
   pars(30) = Yu(0,0);
   pars(31) = Yu(0,1);
   pars(32) = Yu(0,2);
   pars(33) = Yu(1,0);
   pars(34) = Yu(1,1);
   pars(35) = Yu(1,2);
   pars(36) = Yu(2,0);
   pars(37) = Yu(2,1);
   pars(38) = Yu(2,2);
   pars(39) = MuPr;
   pars(40) = MuPhi;
   pars(41) = XiF;
   pars(42) = g1;
   pars(43) = g2;
   pars(44) = g3;
   pars(45) = g1p;
   pars(46) = vd;
   pars(47) = vu;
   pars(48) = vs;
   pars(49) = vsb;
   pars(50) = vphi;


   return pars;
}

void NE6SSM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters:\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "KappaPr = " << KappaPr << '\n';
   ostr << "Sigmax = " << Sigmax << '\n';
   ostr << "Kappa = " << Kappa << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "MuPr = " << MuPr << '\n';
   ostr << "MuPhi = " << MuPhi << '\n';
   ostr << "XiF = " << XiF << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "g1p = " << g1p << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
   ostr << "vs = " << vs << '\n';
   ostr << "vsb = " << vsb << '\n';
   ostr << "vphi = " << vphi << '\n';

}

void NE6SSM_susy_parameters::set(const Eigen::ArrayXd& v)
{
   Yd(0,0) = v(0);
   Yd(0,1) = v(1);
   Yd(0,2) = v(2);
   Yd(1,0) = v(3);
   Yd(1,1) = v(4);
   Yd(1,2) = v(5);
   Yd(2,0) = v(6);
   Yd(2,1) = v(7);
   Yd(2,2) = v(8);
   Ye(0,0) = v(9);
   Ye(0,1) = v(10);
   Ye(0,2) = v(11);
   Ye(1,0) = v(12);
   Ye(1,1) = v(13);
   Ye(1,2) = v(14);
   Ye(2,0) = v(15);
   Ye(2,1) = v(16);
   Ye(2,2) = v(17);
   KappaPr = v(18);
   Sigmax = v(19);
   Kappa(0,0) = v(20);
   Kappa(0,1) = v(21);
   Kappa(0,2) = v(22);
   Kappa(1,0) = v(23);
   Kappa(1,1) = v(24);
   Kappa(1,2) = v(25);
   Kappa(2,0) = v(26);
   Kappa(2,1) = v(27);
   Kappa(2,2) = v(28);
   Lambdax = v(29);
   Yu(0,0) = v(30);
   Yu(0,1) = v(31);
   Yu(0,2) = v(32);
   Yu(1,0) = v(33);
   Yu(1,1) = v(34);
   Yu(1,2) = v(35);
   Yu(2,0) = v(36);
   Yu(2,1) = v(37);
   Yu(2,2) = v(38);
   MuPr = v(39);
   MuPhi = v(40);
   XiF = v(41);
   g1 = v(42);
   g2 = v(43);
   g3 = v(44);
   g1p = v(45);
   vd = v(46);
   vu = v(47);
   vs = v(48);
   vsb = v(49);
   vphi = v(50);

}

const NE6SSM_input_parameters& NE6SSM_susy_parameters::get_input() const
{
   return input;
}

void NE6SSM_susy_parameters::set_input_parameters(const NE6SSM_input_parameters& input_)
{
   input = input_;
}

void NE6SSM_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappa = (Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = (Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;

}

std::ostream& operator<<(std::ostream& ostr, const NE6SSM_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
