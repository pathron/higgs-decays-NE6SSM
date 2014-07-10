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

void higgs_decay_example(NE6SSM<Two_scale>& ne6ssm, bool speak) {
   //print to see what we have
   if(speak) ne6ssm.print(std::cout);

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
   Eigen::Array<double,5,1> mhiggs_pole;
   Eigen::Array<double,5,1> mAhiggs_pole;

   int A = 2;  //index for lightest physical pseudo scalar
   //(0 and 1 are goldstones)
   mhiggs = ne6ssm.get_Mhh();
   mAhiggs = ne6ssm.get_MAh();
   Uhiggs = ne6ssm.get_ZH();
   MHmatrix = ne6ssm.get_mass_matrix_hh();
   //   if(fabs(mAhiggs(0) -  90.882)/ mAhiggs(0)  < 1e-3) A = 1;
   mhiggs_pole = ne6ssm.get_physical().Mhh;
   mAhiggs_pole = ne6ssm.get_physical().MAh;
   if(speak){
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
   }
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
 if(speak){
 std::cout << "GHA1A1 = " << GHA1A1 << std::endl;
 std::cout << "GhA1A1 = " << (Uhiggs * GHA1A1).transpose() << std::endl;
 std::cout << "GhphysA1A1 =" << GhphysA1A1  << std::endl;
 std::cout << "0.5 * GhphysA1A1 = "  << 0.5 * GhphysA1A1  << std::endl;
 }
 //From table 4 of arXiv::1201.2671 for mh = 126 GeV
 double GamSM = 4.085e-03; 

// Note the GHA1A1 is the Feynman rule -- which is a factor two larger
// than the term in the Lagrangian which Roman calculates
// and coupling \xi_{hAA} which is written in Eq C28a and C28b on p530 
// of Baer & Tata weak scale supersymmetry 
 double GamhA1A1 = Gamma_2body_scalar(mhiggs_pole(0), mAhiggs_pole(A), 
                                      Re(GhphysA1A1(0)));
 double GamTot = GamSM + GamhA1A1;
 double BRhA1A1 = GamhA1A1 / GamTot;
 if(speak){
 std::cout << "GamhA1A1 = " << GamhA1A1 << std::endl;
 std::cout << "GamTot = " << GamTot << std::endl;
 std::cout << "BRhA1A1 =" << BRhA1A1 << std::endl;
 } 
bool BenchMark =false;
 if(BenchMark) {
    cout << "BenchMark Point." << endl;
    cout << "\\begin{tabular}{| c || c | c |}"  << std::endl;
    //just output parameters in  long collumn
    std::cout << "\\lambda                    & \t";  
    std::cout <<   ne6ssm.get_Lambdax() << std::endl;
    
    std::cout << "\\sigma                    & \t";  
    std::cout <<   ne6ssm.get_Sigmax() << std::endl;

    std::cout << "\\kappa                    & \t";  
    std::cout <<   ne6ssm.get_KappaPr() << std::endl;
    
    std::cout << "a_{\\lambda} (GeV) \t  & \t";  
    std::cout <<   ne6ssm.get_TLambdax() << std::endl;

    std::cout << "a_{\\sigma} (GeV) \t & \t";  
    std::cout <<   ne6ssm.get_TSigmax() << std::endl;
    
    std::cout << "a_{\\kappa} (GeV) \t & \t";  
    std::cout <<   ne6ssm.get_TKappaPr() << std::endl;

    std::cout << "\\varphi                    & \t";
    std::cout <<   ne6ssm.get_vphi() << std::endl;
    
    std::cout << "s                   \t & \t";
    std::cout <<   Sqrt(Sqr(ne6ssm.get_vsb()) + Sqr(ne6ssm.get_vs())) << std::endl;

    std::cout << "\\tan\\theta                    & \t";
    std::cout <<   ne6ssm.get_vsb() / ne6ssm.get_vs() << std::endl;

    std::cout << "M_1 = M_1^\prime (GeV) \t & \t";
    std::cout <<   ne6ssm.get_MassB() << std::endl;

    std::cout << "m_{\\chi_1^0} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MChi(0) << std::endl;
    
    std::cout << "m_S^2 (GeV)                    & \t";
    std::cout <<   ne6ssm.get_ms2() << std::endl;

     std::cout << "m_\\overline{S}^2 (GeV)                    & \t";
    std::cout <<   ne6ssm.get_msbar2() << std::endl;

     std::cout << "m_\\phi^2 (GeV)                    & \t";
    std::cout <<   ne6ssm.get_mphi2() << std::endl;
    
     std::cout << "m_{H_d}^2 (GeV)                    & \t";
    std::cout <<   ne6ssm.get_mHd2() << std::endl;
    
    std::cout << "m_{H_u}^2 (GeV)                    & \t";
    std::cout <<   ne6ssm.get_mHu2() << std::endl;
    
    std::cout << "m_{Z^\\prime} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MVZp << std::endl;

    std::cout << "m_{H^\\pm} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MHpm(1) << std::endl;

    std::cout << "m_{A_1} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(2) << std::endl;

    std::cout << "m_{A_2} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(3) << std::endl;

    std::cout << "m_{A_3} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(3) << std::endl;
    
    std::cout << "m_{h_1} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(0) << std::endl;

    std::cout << "m_{h_2} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(1) << std::endl;

    std::cout << "m_{h_3} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(2) << std::endl;

     std::cout << "m_{h_4} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(3) << std::endl;

    std::cout << "m_{h_5} (GeV)                    & \t";
    std::cout <<   ne6ssm.get_physical().MAh(4) << std::endl;

    std::cout << "G (GeV)                    & \t";
    std::cout <<   0.5 *  Re(GhphysA1A1(0)) << std::endl;

    std::cout << "BR(h_1\\longarrow A_1A_1)  & \t";
    std::cout << BRhA1A1 << std::endl;

    std::cout << "\\Gamma(h_1\\longarrow A_1A_1) (GeV) & \t";
    std::cout <<  GamhA1A1 << std::endl;

    std::cout << "\\Gamma_{total} (GeV) & \t";
    std::cout <<  GamTot << std::endl;

 }

 cout << "   "  << BRhA1A1 << std::endl;


}


int main(int argc, const char *argv[])
{
   bool speak = false;
   if(speak){
   INFO("=============================");
   INFO("running higgs_decay_example()");
   INFO("=============================");
   }
   GetSpec gs;
   NE6SSM<Two_scale> ne6ssm;
   // NE6SSM<Two_scale> ne6ssm;
   NE6SSM_slha_io slha;

   gs.findSpectrum(argc,argv,ne6ssm, slha);
   
   higgs_decay_example(ne6ssm,speak);
   
   if(speak) std::cout << "End of calculation."  << std::endl;
   
   return 0;
}
