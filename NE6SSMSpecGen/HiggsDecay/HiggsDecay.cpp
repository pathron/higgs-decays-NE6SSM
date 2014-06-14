//This code has been produced using the routines of 
// FlexibleSUSY arXiv:1406.2319 as part of work for arXiv:yymm.xxxx
//Please cite both papers (and the codes which FlexibleSUSY depends on)
//if you use ths code.


#include "NE6SSM_two_scale_model.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

void setup(NE6SSM<Two_scale>& ne6ssm)
{
   Eigen::Matrix<double,3,3> Yu;
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;
   double Mu;
   double g1;
   double g2;
   double g3;
   double vd;
   double vu;
   Eigen::Matrix<double,3,3> TYu;
   Eigen::Matrix<double,3,3> TYd;
   Eigen::Matrix<double,3,3> TYe;
   double BMu;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double MassB;
   double MassWB;
   double MassG;

   // susy parameters
   Yu << 1.26136e-05, 0, 0,
                   0, 0.00667469, 0,
                   0, 0, 0.857849;

   Yd << 0.000242026, 0, 0,
                   0, 0.00529911, 0,
                   0, 0, 0.193602;

   Ye << 2.84161e-05, 0, 0,
                   0, 0.00587557, 0,
                   0, 0, 0.10199;

   Mu = 627.164;
   g1 = 0.468171;
   g2 = 0.642353;
   g3 = 1.06459;
   vd = 25.0944;
   vu = 242.968;

   // soft parameters
   TYu << -0.0144387, 0, 0,
                   0, -7.64037, 0,
                   0, 0, -759.305;

   TYd << -0.336207, 0, 0,
                  0, -7.36109, 0,
                  0, 0, -250.124;

   TYe << -0.00825134, 0, 0,
                    0, -1.70609, 0,
                    0, 0, -29.4466;

   BMu = 52140.8;

   mq2 << 1.03883e+06, 0, 0,
                    0, 1.03881e+06, 0,
                    0, 0, 879135;

   ml2 << 124856, 0, 0,
               0, 124853, 0,
               0, 0, 124142;

   mHd2 = 92436.9;
   mHu2 = -380337;

   md2 << 954454, 0, 0,
               0, 954439, 0,
               0, 0, 934727;

   mu2 << 963422, 0, 0,
               0, 963400, 0,
               0, 0, 656621;

   me2 << 49215.8, 0, 0,
                0, 49210.9, 0,
                0, 0, 47759.2;

   MassB = 210.328;
   MassWB = 389.189;
   MassG = 1114.45;

   // set parameters
   ne6ssm.set_scale(Electroweak_constants::MZ);
   ne6ssm.set_Yu(Yu);
   ne6ssm.set_Yd(Yd);
   ne6ssm.set_Ye(Ye);
   //ne6ssm.set_Mu(Mu);
   ne6ssm.set_g1(g1);
   ne6ssm.set_g2(g2);
   ne6ssm.set_g3(g3);
   ne6ssm.set_vd(vd);
   ne6ssm.set_vu(vu);
   ne6ssm.set_TYu(TYu);
   ne6ssm.set_TYd(TYd);
   ne6ssm.set_TYe(TYe);
   //ne6ssm.set_BMu(BMu);
   ne6ssm.set_mq2(mq2);
   ne6ssm.set_ml2(ml2);
   ne6ssm.set_mHd2(mHd2);
   ne6ssm.set_mHu2(mHu2);
   ne6ssm.set_md2(md2);
   ne6ssm.set_mu2(mu2);
   ne6ssm.set_me2(me2);
   ne6ssm.set_MassB(MassB);
   ne6ssm.set_MassWB(MassWB);
   ne6ssm.set_MassG(MassG);
}


void higgs_decay_example() {

   NE6SSM<Two_scale> ne6ssm;
   setup(ne6ssm);
   ne6ssm.calculate_DRbar_parameters();



}


int main()
{
   INFO("=============================");
   INFO("running higgs_decay_example()");
   INFO("=============================");

   higgs_decay_example();
   return 0;
}
