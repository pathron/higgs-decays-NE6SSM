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
   TYu <<  -1.78213613E-02, 0, 0,
                   0, -9.43028819E+00 , 0,
                   0, 0, -8.84140948E+02 ;

   TYd << -4.34942067E-01, 0, 0,
                  0, -9.52281881E+00, 0,
                  0, 0, -3.17568217E+02;

   TYe << -8.72744876E-03, 0, 0,
                    0, -1.80452101E+00 , 0,
                    0, 0, -3.11090181E+01;


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

   double xiF = 0.0;
   double muPhi = 0.0;
   double kapPr = 0.03;
   double sigPr = 0.1;
   Eigen::Matrix<double,3,3> kap;

   kap << 0.5,0,   0,
          0,  0.5, 0,
          0,  0,   0.5;  
   
   Eigen::Matrix<double,3,3> Tkap;

   Tkap << 500,0,   0,
          0,  500, 0,
      0,  0,   500; 

   Eigen::Matrix<double,3,3>  mDx2;

   mDx2 << 1.29059458E+06, 0.00000000E+00, 0.00000000E+00,   
      0.00000000E+00,  1.29056775E+06, 0.00000000E+00,   
      0.00000000E+00,  0.00000000E+00, 1.08191301E+06;

    Eigen::Matrix<double,3,3>  mDxBar2;

   mDxBar2 << 1.29059458E+06, 0.00000000E+00, 0.00000000E+00,   
      0.00000000E+00,  1.29056775E+06, 0.00000000E+00,   
      0.00000000E+00,  0.00000000E+00, 1.08191301E+06;

   

      

   double lam =0.1;
   double vs = 5916.8425;
   double vsb = 5384.32672;
   double vphi = 6000;
   double gN = 0.47;
   double MuPr = 10000.0;
   double LXiF = 0.0;
   // double ms2 This will be fixed by EWSB  

   double TkapPr = 27.9;
   double Tsig = 120.0;
   double Tlam = 60.0;
   double BMuPr = 10000.0;
   double BMuPhi = 0.0;
   double mHp2 = 1.0E+08;
   double mHpBar2 = 1.0E+08;
   double MassBp = 600.0;

   ne6ssm.set_XiF(xiF);
   ne6ssm.set_MuPhi(muPhi);
   ne6ssm.set_KappaPr(kapPr);
   ne6ssm.set_Sigmax(sigPr);
   ne6ssm.set_Kappa(kap);
   ne6ssm.set_Lambdax(lam);
   ne6ssm.set_vs(vs);
   ne6ssm.set_vsb(vsb);
   ne6ssm.set_vphi(vphi);
   ne6ssm.set_g1p(gN);
   ne6ssm.set_MuPr(MuPr);

   ne6ssm.set_TKappaPr(TkapPr);
   ne6ssm.set_TSigmax(Tsig);
   ne6ssm.set_TKappa(Tkap);
   ne6ssm.set_TLambdax(Tlam);
   ne6ssm.set_BMuPr(BMuPr);
   ne6ssm.set_BMuPhi(BMuPhi);
   ne6ssm.set_LXiF(LXiF);
   //ne6ssm.set_ms2(ms2); This will be fixed by EWSB  
   ne6ssm.set_mDx2(mDx2);
   ne6ssm.set_mDxbar2(mDxBar2);
   ne6ssm.set_mHp2(mHp2); 
   ne6ssm.set_mHpbar2(mHpBar2); 
   ne6ssm.set_MassBp(MassBp);


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
   NE6SSM_input_parameters input;
   input.QS = 5;
   NE6SSM<Two_scale> ne6ssm(input);
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
