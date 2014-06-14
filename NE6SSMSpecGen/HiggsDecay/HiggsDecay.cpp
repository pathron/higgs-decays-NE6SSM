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
   Yu <<  1.25594045E-05, 0, 0,
                   0, 6.64603602E-03, 0,
                   0, 0, 9.01260765E-01;

   Yd << 2.39287126E-04, 0, 0,
                   0, 5.23914003E-03, 0,
                   0, 0, 2.14811963E-01;

   Ye << 2.88373643E-05, 0, 0,
                   0, 5.96266175E-03, 0,
                   0, 0, 1.00328718E-01;

   g1 = 3.70201899E-01/sqrt(0.6);
   g2 = 6.48596334E-01;
   g3 = 1.23250696E+00;
   vd = 2.46858401E+01;
   vu = 2.37406676E+02;

   std::cout << "sqrt(vu*vu + vd*vd) = " <<  sqrt(vu*vu + vd*vd) <<std::endl;
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


   mq2 << 1.29059458E+06, 0, 0,
                    0, 1.29056775E+06  , 0,
                    0, 0, 1.08191301E+06;

   ml2 << 1.29705183E+05, 0, 0,
               0, 1.29705183E+05, 0,
               0, 0, 1.28885261E+05;

 

   md2 << 1.20425905E+06, 0, 0,
               0, 1.20425905E+06, 0,
               0, 0, 1.17488833E+06;

   mu2 << 1.21304977E+06, 0, 0,
               0, 1.21304977E+06, 0,
               0, 0, 8.15783420E+05;

   me2 << 4.97846223E+04, 0, 0,
      0, 4.97846223E+04, 0,
      0, 0, 4.81130118E+04;

   MassB = 600.0;
   MassWB = 1200.0;
   MassG = 1.22669740E+03;

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
   // ne6ssm.set_mHd2(mHd2); set by EWSB
   // ne6ssm.set_mHu2(mHu2); set by EWSB
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
   // apply EWSB constraint
   ne6ssm.solve_ewsb_tree_level();
   // calculate tree-level spectrum
   ne6ssm.calculate_DRbar_parameters();
   //Print tree spectrum for comparison
   ne6ssm.print(std::cout);
   // get hAA coupling
   // Need to know:
   // mass ordering to avoid the goldstone boson 
   // this will vary as aim is to lower the cp odd mass such
   // that light higgs will decay into pseudo scalars.
   //   Eigen::Array<double,3,1> Ggauge;
   //Eigen::VectorXd Ggauge(3);
   
   Eigen::Array<double,5,1> mhiggs;
   Eigen::Array<double,5,1> mAhiggs;
   Eigen::Matrix<double,5,5> Uhiggs;
   Eigen::Matrix<double,5,5> MHmatrix;
   
   mhiggs = ne6ssm.get_Mhh();
   mAhiggs = ne6ssm.get_MAh();
   Uhiggs =  ne6ssm.get_ZH();
   MHmatrix =  ne6ssm.get_mass_matrix_hh();
   
   std::cout << "mhiggs = " << mhiggs << std::endl;
   std::cout << "mAhiggs = " << mAhiggs << std::endl;
   
   
   std::cout << "Diagonalise 1: " << std::endl;
   std::cout << Uhiggs.transpose() *  MHmatrix * Uhiggs
             << std::endl;
   std::cout << "Diagonalise 2: "  << std::endl;
   std::cout << Uhiggs *  MHmatrix * Uhiggs.transpose()
             << std::endl;
   
   
   for(int j=0; j<=4; j++){
      for(int i=0; i<=4; i++){
         std::cout << "ne6ssm.CpUhhAhAh(" <<i << "," <<j << "," << j << ") = " << ne6ssm.CpUhhAhAh(i,j,j) <<std::endl;
      }
   }
   // std::cout << "Ggauge = " << Ggauge << std::endl;
   


// For *this* point we want Ah(1) not Ah(0).
// From above Diagonalisatio2 applies so the convention is h = UH
// where h is mass egenstate and H is gauge eigenstate
// Can write gauge eigenstate couplings with j as G^H dot h A_jA_j(no sum on j)
//                                                = (G^H)^T h A_j A_j(no j sum)
// So (G^H)^T HAA = (G^H)^T U^T hAA = (U G^H)^T hAA = G hAA
Eigen::Matrix<std::complex<double>,5,1> GHA1A1;
for(int j=0; j<=4; j++){
   GHA1A1(j) = ne6ssm.CpUhhAhAh(j,1,1);
   }

 std::cout << "GHA1A1 = " << GHA1A1 << std::endl;
 std::cout << "GhA1A1 = " << (Uhiggs * GHA1A1).transpose() << std::endl;
 

}

int main()
{
   INFO("=============================");
   INFO("running higgs_decay_example()");
   INFO("=============================");

   higgs_decay_example();
   return 0;
}
