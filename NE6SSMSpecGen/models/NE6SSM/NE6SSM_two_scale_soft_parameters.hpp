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

#ifndef NE6SSM_TWO_SCALE_soft_parameters_H
#define NE6SSM_TWO_SCALE_soft_parameters_H

#include "rge.h"
#include "NE6SSM_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class NE6SSM_soft_parameters : public NE6SSM_susy_parameters {
public:
   explicit NE6SSM_soft_parameters(const NE6SSM_input_parameters& input_ = NE6SSM_input_parameters());
   NE6SSM_soft_parameters(const NE6SSM_susy_parameters& , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, double TKappaPr_, double TSigmax_, const Eigen::Matrix<double,3,3>&
   TKappa_, double TLambdax_, const Eigen::Matrix<double,3,3>& TYu_, double
   BMuPr_, double BMuPhi_, double LXiF_, const Eigen::Matrix<double,3,3>& mq2_,
   const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_, const
   Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_,
   const Eigen::Matrix<double,3,3>& me2_, double ms2_, double msbar2_, const
   Eigen::Matrix<double,3,3>& mDx2_, const Eigen::Matrix<double,3,3>& mDxbar2_,
   double mHp2_, double mHpbar2_, double mphi2_, double MassB_, double MassWB_
   , double MassG_, double MassBp_
);
   virtual ~NE6SSM_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual const Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   NE6SSM_soft_parameters calc_beta() const;
   virtual void clear();

   void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) { TYd = TYd_; }
   void set_TYd(int i, int k, double value) { TYd(i,k) = value; }
   void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) { TYe = TYe_; }
   void set_TYe(int i, int k, double value) { TYe(i,k) = value; }
   void set_TKappaPr(double TKappaPr_) { TKappaPr = TKappaPr_; }
   void set_TSigmax(double TSigmax_) { TSigmax = TSigmax_; }
   void set_TKappa(const Eigen::Matrix<double,3,3>& TKappa_) { TKappa = TKappa_; }
   void set_TKappa(int i, int k, double value) { TKappa(i,k) = value; }
   void set_TLambdax(double TLambdax_) { TLambdax = TLambdax_; }
   void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) { TYu = TYu_; }
   void set_TYu(int i, int k, double value) { TYu(i,k) = value; }
   void set_BMuPr(double BMuPr_) { BMuPr = BMuPr_; }
   void set_BMuPhi(double BMuPhi_) { BMuPhi = BMuPhi_; }
   void set_LXiF(double LXiF_) { LXiF = LXiF_; }
   void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) { mq2 = mq2_; }
   void set_mq2(int i, int k, double value) { mq2(i,k) = value; }
   void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) { ml2 = ml2_; }
   void set_ml2(int i, int k, double value) { ml2(i,k) = value; }
   void set_mHd2(double mHd2_) { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) { mHu2 = mHu2_; }
   void set_md2(const Eigen::Matrix<double,3,3>& md2_) { md2 = md2_; }
   void set_md2(int i, int k, double value) { md2(i,k) = value; }
   void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) { mu2 = mu2_; }
   void set_mu2(int i, int k, double value) { mu2(i,k) = value; }
   void set_me2(const Eigen::Matrix<double,3,3>& me2_) { me2 = me2_; }
   void set_me2(int i, int k, double value) { me2(i,k) = value; }
   void set_ms2(double ms2_) { ms2 = ms2_; }
   void set_msbar2(double msbar2_) { msbar2 = msbar2_; }
   void set_mDx2(const Eigen::Matrix<double,3,3>& mDx2_) { mDx2 = mDx2_; }
   void set_mDx2(int i, int k, double value) { mDx2(i,k) = value; }
   void set_mDxbar2(const Eigen::Matrix<double,3,3>& mDxbar2_) { mDxbar2 = mDxbar2_; }
   void set_mDxbar2(int i, int k, double value) { mDxbar2(i,k) = value; }
   void set_mHp2(double mHp2_) { mHp2 = mHp2_; }
   void set_mHpbar2(double mHpbar2_) { mHpbar2 = mHpbar2_; }
   void set_mphi2(double mphi2_) { mphi2 = mphi2_; }
   void set_MassB(double MassB_) { MassB = MassB_; }
   void set_MassWB(double MassWB_) { MassWB = MassWB_; }
   void set_MassG(double MassG_) { MassG = MassG_; }
   void set_MassBp(double MassBp_) { MassBp = MassBp_; }

   const Eigen::Matrix<double,3,3>& get_TYd() const { return TYd; }
   double get_TYd(int i, int k) const { return TYd(i,k); }
   const Eigen::Matrix<double,3,3>& get_TYe() const { return TYe; }
   double get_TYe(int i, int k) const { return TYe(i,k); }
   double get_TKappaPr() const { return TKappaPr; }
   double get_TSigmax() const { return TSigmax; }
   const Eigen::Matrix<double,3,3>& get_TKappa() const { return TKappa; }
   double get_TKappa(int i, int k) const { return TKappa(i,k); }
   double get_TLambdax() const { return TLambdax; }
   const Eigen::Matrix<double,3,3>& get_TYu() const { return TYu; }
   double get_TYu(int i, int k) const { return TYu(i,k); }
   double get_BMuPr() const { return BMuPr; }
   double get_BMuPhi() const { return BMuPhi; }
   double get_LXiF() const { return LXiF; }
   const Eigen::Matrix<double,3,3>& get_mq2() const { return mq2; }
   double get_mq2(int i, int k) const { return mq2(i,k); }
   const Eigen::Matrix<double,3,3>& get_ml2() const { return ml2; }
   double get_ml2(int i, int k) const { return ml2(i,k); }
   double get_mHd2() const { return mHd2; }
   double get_mHu2() const { return mHu2; }
   const Eigen::Matrix<double,3,3>& get_md2() const { return md2; }
   double get_md2(int i, int k) const { return md2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mu2() const { return mu2; }
   double get_mu2(int i, int k) const { return mu2(i,k); }
   const Eigen::Matrix<double,3,3>& get_me2() const { return me2; }
   double get_me2(int i, int k) const { return me2(i,k); }
   double get_ms2() const { return ms2; }
   double get_msbar2() const { return msbar2; }
   const Eigen::Matrix<double,3,3>& get_mDx2() const { return mDx2; }
   double get_mDx2(int i, int k) const { return mDx2(i,k); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2() const { return mDxbar2; }
   double get_mDxbar2(int i, int k) const { return mDxbar2(i,k); }
   double get_mHp2() const { return mHp2; }
   double get_mHpbar2() const { return mHpbar2; }
   double get_mphi2() const { return mphi2; }
   double get_MassB() const { return MassB; }
   double get_MassWB() const { return MassWB; }
   double get_MassG() const { return MassG; }
   double get_MassBp() const { return MassBp; }


protected:
   Eigen::Matrix<double,3,3> TYd;
   Eigen::Matrix<double,3,3> TYe;
   double TKappaPr;
   double TSigmax;
   Eigen::Matrix<double,3,3> TKappa;
   double TLambdax;
   Eigen::Matrix<double,3,3> TYu;
   double BMuPr;
   double BMuPhi;
   double LXiF;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double ms2;
   double msbar2;
   Eigen::Matrix<double,3,3> mDx2;
   Eigen::Matrix<double,3,3> mDxbar2;
   double mHp2;
   double mHpbar2;
   double mphi2;
   double MassB;
   double MassWB;
   double MassG;
   double MassBp;


private:
   static const int numberOfParameters = 167;

   struct Soft_traces {
      double traceAdjYdTYd;
      double traceAdjYeTYe;
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceKappaAdjKappa;
      double traceAdjYuTYu;
      double traceAdjKappaTKappa;
      double traceYdAdjYdTYdAdjYd;
      double traceYdAdjYuTYuAdjYd;
      double traceYeAdjYeTYeAdjYe;
      double traceYuAdjYdTYdAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceKappaAdjKappaKappaAdjKappa;
      double traceKappaAdjKappaTKappaAdjKappa;
      double traceYuAdjYuYuAdjYu;
      double traceYuAdjYuTYuAdjYu;
      double traceconjTYdTpTYd;
      double traceconjTYeTpTYe;
      double tracemd2YdAdjYd;
      double traceme2YeAdjYe;
      double traceml2AdjYeYe;
      double tracemq2AdjYdYd;
      double traceconjTYdTpYd;
      double traceconjTYeTpYe;
      double traceconjTYuTpTYu;
      double tracemq2AdjYuYu;
      double tracemu2YuAdjYu;
      double traceconjTYuTpYu;
      double traceconjTKappaTpKappa;
      double traceconjTKappaTpTKappa;
      double traceKappaAdjKappaconjmDx2;
      double traceKappaconjmDxbar2AdjKappa;
      double traceYdAdjYdTYdAdjTYd;
      double traceYdAdjYuTYuAdjTYd;
      double traceYdAdjTYdTYdAdjYd;
      double traceYdAdjTYuTYuAdjYd;
      double traceYeAdjYeTYeAdjTYe;
      double traceYeAdjTYeTYeAdjYe;
      double traceYuAdjYdTYdAdjTYu;
      double traceYuAdjTYdTYdAdjYu;
      double tracemd2YdAdjYdYdAdjYd;
      double tracemd2YdAdjYuYuAdjYd;
      double traceme2YeAdjYeYeAdjYe;
      double traceml2AdjYeYeAdjYeYe;
      double tracemq2AdjYdYdAdjYdYd;
      double tracemq2AdjYdYdAdjYuYu;
      double tracemq2AdjYuYuAdjYdYd;
      double tracemu2YuAdjYdYdAdjYu;
      double traceYuAdjYuTYuAdjTYu;
      double traceYuAdjTYuTYuAdjYu;
      double tracemq2AdjYuYuAdjYuYu;
      double tracemu2YuAdjYuYuAdjYu;
      double traceKappaAdjKappaTKappaAdjTKappa;
      double traceKappaAdjTKappaTKappaAdjKappa;
      double traceKappaAdjKappaKappaAdjKappaconjmDx2;
      double traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
      double traceKappaAdjKappaconjmDx2KappaAdjKappa;
      double traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
      double Tr11;
      double Tr14;
      double Tr2U111;
      double Tr2U114;
      double Tr31;
      double Tr22;
      double Tr23;
      double Tr2U141;
      double Tr2U144;
      double Tr34;

   };
   void calc_soft_traces(Soft_traces&) const;

   Eigen::Matrix<double,3,3> calc_beta_TYd_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYd_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYe_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappaPr_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TKappaPr_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TSigmax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TSigmax_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TKappa_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TKappa_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_TLambdax_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_TYu_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPr_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPr_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPhi_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_BMuPhi_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_LXiF_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_LXiF_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mq2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_ml2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHd2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_md2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_me2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_ms2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_msbar2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_msbar2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDx2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDx2_two_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDxbar2_one_loop(const TRACE_STRUCT_TYPE&) const;
   Eigen::Matrix<double,3,3> calc_beta_mDxbar2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHp2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHp2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHpbar2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mHpbar2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mphi2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mphi2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassWB_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassG_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassBp_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_MassBp_two_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const NE6SSM_soft_parameters&);

} // namespace flexiblesusy

#endif
