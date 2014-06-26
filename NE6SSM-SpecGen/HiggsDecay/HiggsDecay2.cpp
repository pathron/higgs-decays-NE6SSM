//This code has been produced using the routines of 
// FlexibleSUSY arXiv:1406.2319 as part of work for arXiv:yymm.xxxx
//Please cite both papers (and the codes which FlexibleSUSY depends on)
//if you use ths code.


#include "NE6SSM_two_scale_model.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "wrappers.hpp"
#include "GetSpec.hpp"
#include "NE6SSM_slha_io.hpp"

using namespace flexiblesusy;

double Gamma_2body_scalar(double m1, double m2, double coup) {
   if (m2 > 0.5 * m1) return 0;

   double rt = sqrt(1 - 4.0 *  m2 * m2 / ( m1 * m1));
   double pf = coup * coup / (32.0 * PI * m1);
   return rt*pf;
}

void higgs_decay_example(NE6SSM<Two_scale>& ne6ssm) {
   //print to see what we have
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
   
   int A = 2;  //index for lightest physical pseudo scalar
   //(0 and 1 are goldstones)
   mhiggs = ne6ssm.get_Mhh();
   mAhiggs = ne6ssm.get_MAh();
   Uhiggs = ne6ssm.get_ZH();
   MHmatrix = ne6ssm.get_mass_matrix_hh();
   //   if(fabs(mAhiggs(0) -  90.882)/ mAhiggs(0)  < 1e-3) A = 1;
   
   std::cout << "mhiggs = " << mhiggs << std::endl;
   std::cout << "mAhiggs = " << mAhiggs << std::endl;
   std::cout << "Uhiggs = " << Uhiggs << std::endl;
   
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
 
   std::cout << "A = " <<A << std::endl;  
// For *this* point we want Ah(1) not Ah(0).
// From above Diagonalisatio2 applies so the convention is h = UH
// where h is mass egenstate and H is gauge eigenstate
// Can write gauge eigenstate couplings with j as G^H dot h A_jA_j(no sum on j)
//                                                = (G^H)^T h A_j A_j(no j sum)
// So (G^H)^T HAA = (G^H)^T U^T hAA = (U G^H)^T hAA = G hAA
   Eigen::Matrix<std::complex<double>,5,1> GHA1A1; //gauge basis
   Eigen::Matrix<std::complex<double>,5,1> GhphysA1A1; //mass basis
for(int j=0; j<=4; j++){
   GHA1A1(j) = ne6ssm.CpUhhAhAh(j,A,A);
 }
 GhphysA1A1 = Uhiggs * GHA1A1;
 std::cout << "GHA1A1 = " << GHA1A1 << std::endl;
 std::cout << "GhA1A1 = " << (Uhiggs * GHA1A1).transpose() << std::endl;
 std::cout << "GhphysA1A1 =" << GhphysA1A1  << std::endl;
 std::cout << "0.5 * GhphysA1A1 = "  << 0.5 * GhphysA1A1  << std::endl;
 //From table 4 of arXiv::1201.2671 for mh = 126 GeV
 double GamSM = 4.085e-03; 

// Note the GHA1A1 is the Feynman rule -- which is a factor two larger
// than the term in the Lagrangian which Roman calculates
// and coupling \xi_{hAA} which is written in Eq C28a and C28b on p530 
// of Baer & Tata weak scale supersymmetry 
 double GamhA1A1 = Gamma_2body_scalar( mhiggs(0), mAhiggs(A), Re(GhphysA1A1(0)));
 double GamTot = GamSM + GamhA1A1;
 double BRhA1A1 = GamhA1A1 / GamTot;
 std::cout << "GamhA1A1 = " << GamhA1A1 << std::endl;
 std::cout << "GamTot = " << GamTot << std::endl;
 std::cout << "BRhA1A1 =" << BRhA1A1 << std::endl;
}


int main(int argc, const char *argv[])
{
   INFO("=============================");
   INFO("running higgs_decay_example()");
   INFO("=============================");
   GetSpec gs;
   NE6SSM<Two_scale> ne6ssm;
   // NE6SSM<Two_scale> ne6ssm;
   NE6SSM_slha_io slha;

   gs.findSpectrum(argc,argv,ne6ssm, slha);

   higgs_decay_example(ne6ssm);
   
   std::cout << "End of calculation."  << std::endl;
   
   return 0;
}
