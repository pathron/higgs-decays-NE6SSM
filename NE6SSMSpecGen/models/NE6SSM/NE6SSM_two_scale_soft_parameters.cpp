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

// File generated at Fri 13 Jun 2014 15:36:36

#include "NE6SSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

NE6SSM_soft_parameters::NE6SSM_soft_parameters(const NE6SSM_input_parameters& input_)
   : NE6SSM_susy_parameters(input_)
   , TYd(Eigen::Matrix<double,3,3>::Zero()), TYe(Eigen::Matrix<double,3,3>
   ::Zero()), TKappaPr(0), TSigmax(0), TKappa(Eigen::Matrix<double,3,3>::Zero()
   ), TLambdax(0), TYu(Eigen::Matrix<double,3,3>::Zero()), BMuPr(0), BMuPhi(0),
   LXiF(0), mq2(Eigen::Matrix<double,3,3>::Zero()), ml2(Eigen::Matrix<double,3
   ,3>::Zero()), mHd2(0), mHu2(0), md2(Eigen::Matrix<double,3,3>::Zero()), mu2(
   Eigen::Matrix<double,3,3>::Zero()), me2(Eigen::Matrix<double,3,3>::Zero()),
   ms2(0), msbar2(0), mDx2(Eigen::Matrix<double,3,3>::Zero()), mDxbar2(
   Eigen::Matrix<double,3,3>::Zero()), mHp2(0), mHpbar2(0), mphi2(0), MassB(0),
   MassWB(0), MassG(0), MassBp(0)

{
   set_number_of_parameters(numberOfParameters);
}

NE6SSM_soft_parameters::NE6SSM_soft_parameters(
   const NE6SSM_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, double TKappaPr_, double TSigmax_, const Eigen::Matrix<double,3,3>&
   TKappa_, double TLambdax_, const Eigen::Matrix<double,3,3>& TYu_, double
   BMuPr_, double BMuPhi_, double LXiF_, const Eigen::Matrix<double,3,3>& mq2_,
   const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_, const
   Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_,
   const Eigen::Matrix<double,3,3>& me2_, double ms2_, double msbar2_, const
   Eigen::Matrix<double,3,3>& mDx2_, const Eigen::Matrix<double,3,3>& mDxbar2_,
   double mHp2_, double mHpbar2_, double mphi2_, double MassB_, double MassWB_
   , double MassG_, double MassBp_

)
   : NE6SSM_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TKappaPr(TKappaPr_), TSigmax(TSigmax_), TKappa(
   TKappa_), TLambdax(TLambdax_), TYu(TYu_), BMuPr(BMuPr_), BMuPhi(BMuPhi_),
   LXiF(LXiF_), mq2(mq2_), ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_), md2(md2_), mu2(
   mu2_), me2(me2_), ms2(ms2_), msbar2(msbar2_), mDx2(mDx2_), mDxbar2(mDxbar2_)
   , mHp2(mHp2_), mHpbar2(mHpbar2_), mphi2(mphi2_), MassB(MassB_), MassWB(
   MassWB_), MassG(MassG_), MassBp(MassBp_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd NE6SSM_soft_parameters::beta() const
{
   return calc_beta().get();
}

NE6SSM_soft_parameters NE6SSM_soft_parameters::calc_beta() const
{
   Soft_traces soft_traces;
   calc_soft_traces(soft_traces);

   Eigen::Matrix<double,3,3> beta_TYd(calc_beta_TYd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYe(calc_beta_TYe_one_loop(TRACE_STRUCT));
   double beta_TKappaPr(calc_beta_TKappaPr_one_loop(TRACE_STRUCT));
   double beta_TSigmax(calc_beta_TSigmax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TKappa(calc_beta_TKappa_one_loop(TRACE_STRUCT
      ));
   double beta_TLambdax(calc_beta_TLambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYu(calc_beta_TYu_one_loop(TRACE_STRUCT));
   double beta_BMuPr(calc_beta_BMuPr_one_loop(TRACE_STRUCT));
   double beta_BMuPhi(calc_beta_BMuPhi_one_loop(TRACE_STRUCT));
   double beta_LXiF(calc_beta_LXiF_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mq2(calc_beta_mq2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_ml2(calc_beta_ml2_one_loop(TRACE_STRUCT));
   double beta_mHd2(calc_beta_mHd2_one_loop(TRACE_STRUCT));
   double beta_mHu2(calc_beta_mHu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_md2(calc_beta_md2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mu2(calc_beta_mu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_me2(calc_beta_me2_one_loop(TRACE_STRUCT));
   double beta_ms2(calc_beta_ms2_one_loop(TRACE_STRUCT));
   double beta_msbar2(calc_beta_msbar2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mDx2(calc_beta_mDx2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mDxbar2(calc_beta_mDxbar2_one_loop(
      TRACE_STRUCT));
   double beta_mHp2(calc_beta_mHp2_one_loop(TRACE_STRUCT));
   double beta_mHpbar2(calc_beta_mHpbar2_one_loop(TRACE_STRUCT));
   double beta_mphi2(calc_beta_mphi2_one_loop(TRACE_STRUCT));
   double beta_MassB(calc_beta_MassB_one_loop(TRACE_STRUCT));
   double beta_MassWB(calc_beta_MassWB_one_loop(TRACE_STRUCT));
   double beta_MassG(calc_beta_MassG_one_loop(TRACE_STRUCT));
   double beta_MassBp(calc_beta_MassBp_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_TYd += calc_beta_TYd_two_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_two_loop(TRACE_STRUCT);
      beta_TKappaPr += calc_beta_TKappaPr_two_loop(TRACE_STRUCT);
      beta_TSigmax += calc_beta_TSigmax_two_loop(TRACE_STRUCT);
      beta_TKappa += calc_beta_TKappa_two_loop(TRACE_STRUCT);
      beta_TLambdax += calc_beta_TLambdax_two_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_two_loop(TRACE_STRUCT);
      beta_BMuPr += calc_beta_BMuPr_two_loop(TRACE_STRUCT);
      beta_BMuPhi += calc_beta_BMuPhi_two_loop(TRACE_STRUCT);
      beta_LXiF += calc_beta_LXiF_two_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_two_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_two_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_two_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_two_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_two_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_two_loop(TRACE_STRUCT);
      beta_ms2 += calc_beta_ms2_two_loop(TRACE_STRUCT);
      beta_msbar2 += calc_beta_msbar2_two_loop(TRACE_STRUCT);
      beta_mDx2 += calc_beta_mDx2_two_loop(TRACE_STRUCT);
      beta_mDxbar2 += calc_beta_mDxbar2_two_loop(TRACE_STRUCT);
      beta_mHp2 += calc_beta_mHp2_two_loop(TRACE_STRUCT);
      beta_mHpbar2 += calc_beta_mHpbar2_two_loop(TRACE_STRUCT);
      beta_mphi2 += calc_beta_mphi2_two_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);
      beta_MassBp += calc_beta_MassBp_two_loop(TRACE_STRUCT);

   }


   const NE6SSM_susy_parameters susy_betas(NE6SSM_susy_parameters::calc_beta());

   return NE6SSM_soft_parameters(susy_betas, beta_TYd, beta_TYe, beta_TKappaPr, beta_TSigmax, beta_TKappa, beta_TLambdax, beta_TYu, beta_BMuPr, beta_BMuPhi, beta_LXiF, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_ms2, beta_msbar2, beta_mDx2, beta_mDxbar2, beta_mHp2, beta_mHpbar2, beta_mphi2, beta_MassB, beta_MassWB, beta_MassG, beta_MassBp);
}

void NE6SSM_soft_parameters::clear()
{
   NE6SSM_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TKappaPr = 0.;
   TSigmax = 0.;
   TKappa = Eigen::Matrix<double,3,3>::Zero();
   TLambdax = 0.;
   TYu = Eigen::Matrix<double,3,3>::Zero();
   BMuPr = 0.;
   BMuPhi = 0.;
   LXiF = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   ms2 = 0.;
   msbar2 = 0.;
   mDx2 = Eigen::Matrix<double,3,3>::Zero();
   mDxbar2 = Eigen::Matrix<double,3,3>::Zero();
   mHp2 = 0.;
   mHpbar2 = 0.;
   mphi2 = 0.;
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;
   MassBp = 0.;

}

const Eigen::ArrayXd NE6SSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(NE6SSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(51) = TYd(0,0);
   pars(52) = TYd(0,1);
   pars(53) = TYd(0,2);
   pars(54) = TYd(1,0);
   pars(55) = TYd(1,1);
   pars(56) = TYd(1,2);
   pars(57) = TYd(2,0);
   pars(58) = TYd(2,1);
   pars(59) = TYd(2,2);
   pars(60) = TYe(0,0);
   pars(61) = TYe(0,1);
   pars(62) = TYe(0,2);
   pars(63) = TYe(1,0);
   pars(64) = TYe(1,1);
   pars(65) = TYe(1,2);
   pars(66) = TYe(2,0);
   pars(67) = TYe(2,1);
   pars(68) = TYe(2,2);
   pars(69) = TKappaPr;
   pars(70) = TSigmax;
   pars(71) = TKappa(0,0);
   pars(72) = TKappa(0,1);
   pars(73) = TKappa(0,2);
   pars(74) = TKappa(1,0);
   pars(75) = TKappa(1,1);
   pars(76) = TKappa(1,2);
   pars(77) = TKappa(2,0);
   pars(78) = TKappa(2,1);
   pars(79) = TKappa(2,2);
   pars(80) = TLambdax;
   pars(81) = TYu(0,0);
   pars(82) = TYu(0,1);
   pars(83) = TYu(0,2);
   pars(84) = TYu(1,0);
   pars(85) = TYu(1,1);
   pars(86) = TYu(1,2);
   pars(87) = TYu(2,0);
   pars(88) = TYu(2,1);
   pars(89) = TYu(2,2);
   pars(90) = BMuPr;
   pars(91) = BMuPhi;
   pars(92) = LXiF;
   pars(93) = mq2(0,0);
   pars(94) = mq2(0,1);
   pars(95) = mq2(0,2);
   pars(96) = mq2(1,0);
   pars(97) = mq2(1,1);
   pars(98) = mq2(1,2);
   pars(99) = mq2(2,0);
   pars(100) = mq2(2,1);
   pars(101) = mq2(2,2);
   pars(102) = ml2(0,0);
   pars(103) = ml2(0,1);
   pars(104) = ml2(0,2);
   pars(105) = ml2(1,0);
   pars(106) = ml2(1,1);
   pars(107) = ml2(1,2);
   pars(108) = ml2(2,0);
   pars(109) = ml2(2,1);
   pars(110) = ml2(2,2);
   pars(111) = mHd2;
   pars(112) = mHu2;
   pars(113) = md2(0,0);
   pars(114) = md2(0,1);
   pars(115) = md2(0,2);
   pars(116) = md2(1,0);
   pars(117) = md2(1,1);
   pars(118) = md2(1,2);
   pars(119) = md2(2,0);
   pars(120) = md2(2,1);
   pars(121) = md2(2,2);
   pars(122) = mu2(0,0);
   pars(123) = mu2(0,1);
   pars(124) = mu2(0,2);
   pars(125) = mu2(1,0);
   pars(126) = mu2(1,1);
   pars(127) = mu2(1,2);
   pars(128) = mu2(2,0);
   pars(129) = mu2(2,1);
   pars(130) = mu2(2,2);
   pars(131) = me2(0,0);
   pars(132) = me2(0,1);
   pars(133) = me2(0,2);
   pars(134) = me2(1,0);
   pars(135) = me2(1,1);
   pars(136) = me2(1,2);
   pars(137) = me2(2,0);
   pars(138) = me2(2,1);
   pars(139) = me2(2,2);
   pars(140) = ms2;
   pars(141) = msbar2;
   pars(142) = mDx2(0,0);
   pars(143) = mDx2(0,1);
   pars(144) = mDx2(0,2);
   pars(145) = mDx2(1,0);
   pars(146) = mDx2(1,1);
   pars(147) = mDx2(1,2);
   pars(148) = mDx2(2,0);
   pars(149) = mDx2(2,1);
   pars(150) = mDx2(2,2);
   pars(151) = mDxbar2(0,0);
   pars(152) = mDxbar2(0,1);
   pars(153) = mDxbar2(0,2);
   pars(154) = mDxbar2(1,0);
   pars(155) = mDxbar2(1,1);
   pars(156) = mDxbar2(1,2);
   pars(157) = mDxbar2(2,0);
   pars(158) = mDxbar2(2,1);
   pars(159) = mDxbar2(2,2);
   pars(160) = mHp2;
   pars(161) = mHpbar2;
   pars(162) = mphi2;
   pars(163) = MassB;
   pars(164) = MassWB;
   pars(165) = MassG;
   pars(166) = MassBp;


   return pars;
}

void NE6SSM_soft_parameters::print(std::ostream& ostr) const
{
   NE6SSM_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TKappaPr = " << TKappaPr << '\n';
   ostr << "TSigmax = " << TSigmax << '\n';
   ostr << "TKappa = " << TKappa << '\n';
   ostr << "TLambdax = " << TLambdax << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "BMuPr = " << BMuPr << '\n';
   ostr << "BMuPhi = " << BMuPhi << '\n';
   ostr << "LXiF = " << LXiF << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "ms2 = " << ms2 << '\n';
   ostr << "msbar2 = " << msbar2 << '\n';
   ostr << "mDx2 = " << mDx2 << '\n';
   ostr << "mDxbar2 = " << mDxbar2 << '\n';
   ostr << "mHp2 = " << mHp2 << '\n';
   ostr << "mHpbar2 = " << mHpbar2 << '\n';
   ostr << "mphi2 = " << mphi2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';
   ostr << "MassBp = " << MassBp << '\n';

}

void NE6SSM_soft_parameters::set(const Eigen::ArrayXd& v)
{
   NE6SSM_susy_parameters::set(v);

   TYd(0,0) = v(51);
   TYd(0,1) = v(52);
   TYd(0,2) = v(53);
   TYd(1,0) = v(54);
   TYd(1,1) = v(55);
   TYd(1,2) = v(56);
   TYd(2,0) = v(57);
   TYd(2,1) = v(58);
   TYd(2,2) = v(59);
   TYe(0,0) = v(60);
   TYe(0,1) = v(61);
   TYe(0,2) = v(62);
   TYe(1,0) = v(63);
   TYe(1,1) = v(64);
   TYe(1,2) = v(65);
   TYe(2,0) = v(66);
   TYe(2,1) = v(67);
   TYe(2,2) = v(68);
   TKappaPr = v(69);
   TSigmax = v(70);
   TKappa(0,0) = v(71);
   TKappa(0,1) = v(72);
   TKappa(0,2) = v(73);
   TKappa(1,0) = v(74);
   TKappa(1,1) = v(75);
   TKappa(1,2) = v(76);
   TKappa(2,0) = v(77);
   TKappa(2,1) = v(78);
   TKappa(2,2) = v(79);
   TLambdax = v(80);
   TYu(0,0) = v(81);
   TYu(0,1) = v(82);
   TYu(0,2) = v(83);
   TYu(1,0) = v(84);
   TYu(1,1) = v(85);
   TYu(1,2) = v(86);
   TYu(2,0) = v(87);
   TYu(2,1) = v(88);
   TYu(2,2) = v(89);
   BMuPr = v(90);
   BMuPhi = v(91);
   LXiF = v(92);
   mq2(0,0) = v(93);
   mq2(0,1) = v(94);
   mq2(0,2) = v(95);
   mq2(1,0) = v(96);
   mq2(1,1) = v(97);
   mq2(1,2) = v(98);
   mq2(2,0) = v(99);
   mq2(2,1) = v(100);
   mq2(2,2) = v(101);
   ml2(0,0) = v(102);
   ml2(0,1) = v(103);
   ml2(0,2) = v(104);
   ml2(1,0) = v(105);
   ml2(1,1) = v(106);
   ml2(1,2) = v(107);
   ml2(2,0) = v(108);
   ml2(2,1) = v(109);
   ml2(2,2) = v(110);
   mHd2 = v(111);
   mHu2 = v(112);
   md2(0,0) = v(113);
   md2(0,1) = v(114);
   md2(0,2) = v(115);
   md2(1,0) = v(116);
   md2(1,1) = v(117);
   md2(1,2) = v(118);
   md2(2,0) = v(119);
   md2(2,1) = v(120);
   md2(2,2) = v(121);
   mu2(0,0) = v(122);
   mu2(0,1) = v(123);
   mu2(0,2) = v(124);
   mu2(1,0) = v(125);
   mu2(1,1) = v(126);
   mu2(1,2) = v(127);
   mu2(2,0) = v(128);
   mu2(2,1) = v(129);
   mu2(2,2) = v(130);
   me2(0,0) = v(131);
   me2(0,1) = v(132);
   me2(0,2) = v(133);
   me2(1,0) = v(134);
   me2(1,1) = v(135);
   me2(1,2) = v(136);
   me2(2,0) = v(137);
   me2(2,1) = v(138);
   me2(2,2) = v(139);
   ms2 = v(140);
   msbar2 = v(141);
   mDx2(0,0) = v(142);
   mDx2(0,1) = v(143);
   mDx2(0,2) = v(144);
   mDx2(1,0) = v(145);
   mDx2(1,1) = v(146);
   mDx2(1,2) = v(147);
   mDx2(2,0) = v(148);
   mDx2(2,1) = v(149);
   mDx2(2,2) = v(150);
   mDxbar2(0,0) = v(151);
   mDxbar2(0,1) = v(152);
   mDxbar2(0,2) = v(153);
   mDxbar2(1,0) = v(154);
   mDxbar2(1,1) = v(155);
   mDxbar2(1,2) = v(156);
   mDxbar2(2,0) = v(157);
   mDxbar2(2,1) = v(158);
   mDxbar2(2,2) = v(159);
   mHp2 = v(160);
   mHpbar2 = v(161);
   mphi2 = v(162);
   MassB = v(163);
   MassWB = v(164);
   MassG = v(165);
   MassBp = v(166);

}

void NE6SSM_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   TRACE_STRUCT.traceAdjYdTYd = (Yd.adjoint()*TYd).trace();
   TRACE_STRUCT.traceAdjYeTYe = (Ye.adjoint()*TYe).trace();
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappa = (Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceAdjYuTYu = (Yu.adjoint()*TYu).trace();
   TRACE_STRUCT.traceAdjKappaTKappa = ((Kappa).adjoint()*TKappa).trace();
   TRACE_STRUCT.traceYdAdjYdTYdAdjYd = (Yd*Yd.adjoint()*TYd*Yd.adjoint()).trace
      ();
   TRACE_STRUCT.traceYdAdjYuTYuAdjYd = (Yd*Yu.adjoint()*TYu*Yd.adjoint()).trace
      ();
   TRACE_STRUCT.traceYeAdjYeTYeAdjYe = (Ye*Ye.adjoint()*TYe*Ye.adjoint()).trace
      ();
   TRACE_STRUCT.traceYuAdjYdTYdAdjYu = (Yu*Yd.adjoint()*TYd*Yu.adjoint()).trace
      ();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = (Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa = (Kappa*(Kappa).adjoint()*
      TKappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYuAdjYuTYuAdjYu = (Yu*Yu.adjoint()*TYu*Yu.adjoint()).trace
      ();
   TRACE_STRUCT.traceconjTYdTpTYd = (TYd.conjugate()*(TYd).transpose()).trace()
      ;
   TRACE_STRUCT.traceconjTYeTpTYe = (TYe.conjugate()*(TYe).transpose()).trace()
      ;
   TRACE_STRUCT.tracemd2YdAdjYd = (md2*Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceme2YeAdjYe = (me2*Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceml2AdjYeYe = (ml2*Ye.adjoint()*Ye).trace();
   TRACE_STRUCT.tracemq2AdjYdYd = (mq2*Yd.adjoint()*Yd).trace();
   TRACE_STRUCT.traceconjTYdTpYd = (TYd.conjugate()*Yd.transpose()).trace();
   TRACE_STRUCT.traceconjTYeTpYe = (TYe.conjugate()*Ye.transpose()).trace();
   TRACE_STRUCT.traceconjTYuTpTYu = (TYu.conjugate()*(TYu).transpose()).trace()
      ;
   TRACE_STRUCT.tracemq2AdjYuYu = (mq2*Yu.adjoint()*Yu).trace();
   TRACE_STRUCT.tracemu2YuAdjYu = (mu2*Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceconjTYuTpYu = (TYu.conjugate()*Yu.transpose()).trace();
   TRACE_STRUCT.traceconjTKappaTpKappa = (TKappa.conjugate()*(Kappa).transpose(
      )).trace();
   TRACE_STRUCT.traceconjTKappaTpTKappa = (TKappa.conjugate()*(TKappa)
      .transpose()).trace();
   TRACE_STRUCT.traceKappaAdjKappaconjmDx2 = (Kappa*(Kappa).adjoint()*
      mDx2.conjugate()).trace();
   TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa = (Kappa*mDxbar2.conjugate()*(
      Kappa).adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYdTYdAdjTYd = (Yd*Yd.adjoint()*TYd*(TYd).adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjYuTYuAdjTYd = (Yd*Yu.adjoint()*TYu*(TYd).adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjTYdTYdAdjYd = (Yd*(TYd).adjoint()*TYd*Yd.adjoint())
      .trace();
   TRACE_STRUCT.traceYdAdjTYuTYuAdjYd = (Yd*(TYu).adjoint()*TYu*Yd.adjoint())
      .trace();
   TRACE_STRUCT.traceYeAdjYeTYeAdjTYe = (Ye*Ye.adjoint()*TYe*(TYe).adjoint())
      .trace();
   TRACE_STRUCT.traceYeAdjTYeTYeAdjYe = (Ye*(TYe).adjoint()*TYe*Ye.adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjYdTYdAdjTYu = (Yu*Yd.adjoint()*TYd*(TYu).adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjTYdTYdAdjYu = (Yu*(TYd).adjoint()*TYd*Yu.adjoint())
      .trace();
   TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = (md2*Yd*Yd.adjoint()*Yd*Yd.adjoint())
      .trace();
   TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = (md2*Yd*Yu.adjoint()*Yu*Yd.adjoint())
      .trace();
   TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = (me2*Ye*Ye.adjoint()*Ye*Ye.adjoint())
      .trace();
   TRACE_STRUCT.traceml2AdjYeYeAdjYeYe = (ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye)
      .trace();
   TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd = (mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd)
      .trace();
   TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu = (mq2*Yd.adjoint()*Yd*Yu.adjoint()*Yu)
      .trace();
   TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd = (mq2*Yu.adjoint()*Yu*Yd.adjoint()*Yd)
      .trace();
   TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu = (mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjYuTYuAdjTYu = (Yu*Yu.adjoint()*TYu*(TYu).adjoint())
      .trace();
   TRACE_STRUCT.traceYuAdjTYuTYuAdjYu = (Yu*(TYu).adjoint()*TYu*Yu.adjoint())
      .trace();
   TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu = (mq2*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      .trace();
   TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu = (mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint())
      .trace();
   TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa = (Kappa*(Kappa).adjoint()*
      TKappa*(TKappa).adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa = (Kappa*(TKappa).adjoint()*
      TKappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2 = (Kappa*(Kappa)
      .adjoint()*Kappa*(Kappa).adjoint()*mDx2.conjugate()).trace();
   TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa = (Kappa*(Kappa)
      .adjoint()*Kappa*mDxbar2.conjugate()*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa = (Kappa*(Kappa)
      .adjoint()*mDx2.conjugate()*Kappa*(Kappa).adjoint()).trace();
   TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa = (Kappa*
      mDxbar2.conjugate()*(Kappa).adjoint()*Kappa*(Kappa).adjoint()).trace();

   const auto QS = INPUT(QS);

   TRACE_STRUCT.Tr11 = 0.7745966692414834*g1*(-mHd2 - mHp2 + mHpbar2 + mHu2 + (
      md2).trace() - (mDx2).trace() + (mDxbar2).trace() + (me2).trace() - (ml2)
      .trace() + (mq2).trace() - 2*(mu2).trace());
   TRACE_STRUCT.Tr14 = 0.15811388300841897*g1p*(-6*mHd2 + 4*mHp2 - 4*mHpbar2 -
      4*mHu2 + ms2*QS - msbar2*QS + 6*(md2).trace() - 6*(mDx2).trace() - 9*(
      mDxbar2).trace() + (me2).trace() + 4*(ml2).trace() + 6*(mq2).trace() + 3*(
      mu2).trace());
   TRACE_STRUCT.Tr2U111 = 0.1*Sqr(g1)*(3*mHd2 + 3*mHp2 + 3*mHpbar2 + 3*mHu2 + 2
      *(md2).trace() + 2*(mDx2).trace() + 2*(mDxbar2).trace() + 6*(me2).trace() +
      3*(ml2).trace() + (mq2).trace() + 8*(mu2).trace());
   TRACE_STRUCT.Tr2U114 = 0.1224744871391589*g1*g1p*(3*mHd2 - 2*mHp2 - 2*
      mHpbar2 - 2*mHu2 + 2*(md2).trace() + 2*(mDx2).trace() - 3*(mDxbar2).trace()
      + (me2).trace() - 2*(ml2).trace() + (mq2).trace() - 2*(mu2).trace());
   TRACE_STRUCT.Tr31 = 0.006454972243679028*g1*(60*(mHd2 - mHu2)*AbsSqr(Lambdax
      ) - 18*mHd2*Sqr(g1) - 18*mHp2*Sqr(g1) + 18*mHpbar2*Sqr(g1) + 18*mHu2*Sqr(g1)
      - 27*mHd2*Sqr(g1p) - 12*mHp2*Sqr(g1p) + 12*mHpbar2*Sqr(g1p) + 12*mHu2*Sqr(
      g1p) - 90*mHd2*Sqr(g2) - 90*mHp2*Sqr(g2) + 90*mHpbar2*Sqr(g2) + 90*mHu2*Sqr(
      g2) + 4*(2*Sqr(g1) + 3*Sqr(g1p) + 40*Sqr(g3))*(md2).trace() - 8*Sqr(g1)*(
      mDx2).trace() - 12*Sqr(g1p)*(mDx2).trace() - 160*Sqr(g3)*(mDx2).trace() + 8*
      Sqr(g1)*(mDxbar2).trace() + 27*Sqr(g1p)*(mDxbar2).trace() + 160*Sqr(g3)*(
      mDxbar2).trace() + 72*Sqr(g1)*(me2).trace() + 3*Sqr(g1p)*(me2).trace() - 18*
      Sqr(g1)*(ml2).trace() - 12*Sqr(g1p)*(ml2).trace() - 90*Sqr(g2)*(ml2).trace()
      + 2*Sqr(g1)*(mq2).trace() + 3*Sqr(g1p)*(mq2).trace() + 90*Sqr(g2)*(mq2)
      .trace() + 160*Sqr(g3)*(mq2).trace() - 64*Sqr(g1)*(mu2).trace() - 6*Sqr(g1p)
      *(mu2).trace() - 320*Sqr(g3)*(mu2).trace() + 180*mHd2*(Yd*Yd.adjoint())
      .trace() + 60*mHd2*(Ye*Ye.adjoint()).trace() - 180*mHu2*(Yu*Yu.adjoint())
      .trace() + 60*(mDx2*Kappa*(Kappa).adjoint()).trace() - 60*(mDxbar2*(Kappa)
      .adjoint()*Kappa).trace() - 120*(Yd*Yd.adjoint()*md2.conjugate()).trace() -
      60*(Yd*mq2.conjugate()*Yd.adjoint()).trace() - 120*(Ye*Ye.adjoint()*
      me2.conjugate()).trace() + 60*(Ye*ml2.conjugate()*Ye.adjoint()).trace() +
      240*(Yu*Yu.adjoint()*mu2.conjugate()).trace() - 60*(Yu*mq2.conjugate()*
      Yu.adjoint()).trace());
   TRACE_STRUCT.Tr22 = 0.5*(mHd2 + mHp2 + mHpbar2 + mHu2 + (ml2).trace() + 3*(
      mq2).trace());
   TRACE_STRUCT.Tr23 = 0.5*((md2).trace() + (mDx2).trace() + (mDxbar2).trace()
      + 2*(mq2).trace() + (mu2).trace());
   TRACE_STRUCT.Tr2U141 = 0.1224744871391589*g1*g1p*(3*mHd2 - 2*mHp2 - 2*
      mHpbar2 - 2*mHu2 + 2*(md2).trace() + 2*(mDx2).trace() - 3*(mDxbar2).trace()
      + (me2).trace() - 2*(ml2).trace() + (mq2).trace() - 2*(mu2).trace());
   TRACE_STRUCT.Tr2U144 = 0.025*Sqr(g1p)*(18*mHd2 + 8*mHp2 + 8*mHpbar2 + 8*mHu2
      + ms2*Sqr(QS) + msbar2*Sqr(QS) + 12*(md2).trace() + 12*(mDx2).trace() + 27*
      (mDxbar2).trace() + (me2).trace() + 8*(ml2).trace() + 6*(mq2).trace() + 3*(
      mu2).trace());
   TRACE_STRUCT.Tr34 = 0.003952847075210475*g1p*(40*(3*mHd2 + 2*mHu2 - ms2*QS)*
      AbsSqr(Lambdax) - 20*(ms2 - msbar2)*QS*AbsSqr(Sigmax) - 36*mHd2*Sqr(g1) + 24
      *mHp2*Sqr(g1) - 24*mHpbar2*Sqr(g1) - 24*mHu2*Sqr(g1) - 54*mHd2*Sqr(g1p) + 16
      *mHp2*Sqr(g1p) - 16*mHpbar2*Sqr(g1p) - 16*mHu2*Sqr(g1p) + ms2*Power(QS,3)*
      Sqr(g1p) - msbar2*Power(QS,3)*Sqr(g1p) - 180*mHd2*Sqr(g2) + 120*mHp2*Sqr(g2)
      - 120*mHpbar2*Sqr(g2) - 120*mHu2*Sqr(g2) + 16*Sqr(g1)*(md2).trace() + 24*
      Sqr(g1p)*(md2).trace() + 320*Sqr(g3)*(md2).trace() - 16*Sqr(g1)*(mDx2).trace
      () - 24*Sqr(g1p)*(mDx2).trace() - 320*Sqr(g3)*(mDx2).trace() - 24*Sqr(g1)*(
      mDxbar2).trace() - 81*Sqr(g1p)*(mDxbar2).trace() - 480*Sqr(g3)*(mDxbar2)
      .trace() + 24*Sqr(g1)*(me2).trace() + Sqr(g1p)*(me2).trace() + 24*Sqr(g1)*(
      ml2).trace() + 16*Sqr(g1p)*(ml2).trace() + 120*Sqr(g2)*(ml2).trace() + 4*Sqr
      (g1)*(mq2).trace() + 6*Sqr(g1p)*(mq2).trace() + 180*Sqr(g2)*(mq2).trace() +
      320*Sqr(g3)*(mq2).trace() + 32*Sqr(g1)*(mu2).trace() + 3*Sqr(g1p)*(mu2)
      .trace() + 160*Sqr(g3)*(mu2).trace() + 360*mHd2*(Yd*Yd.adjoint()).trace() +
      120*mHd2*(Ye*Ye.adjoint()).trace() + 240*mHu2*(Yu*Yu.adjoint()).trace() - 60
      *ms2*QS*(Kappa*(Kappa).adjoint()).trace() + 120*(mDx2*Kappa*(Kappa).adjoint(
      )).trace() + 180*(mDxbar2*(Kappa).adjoint()*Kappa).trace() - 240*(Yd*
      Yd.adjoint()*md2.conjugate()).trace() - 120*(Yd*mq2.conjugate()*Yd.adjoint()
      ).trace() - 40*(Ye*Ye.adjoint()*me2.conjugate()).trace() - 80*(Ye*
      ml2.conjugate()*Ye.adjoint()).trace() - 120*(Yu*Yu.adjoint()*mu2.conjugate()
      ).trace() - 120*(Yu*mq2.conjugate()*Yu.adjoint()).trace());

}

std::ostream& operator<<(std::ostream& ostr, const NE6SSM_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
