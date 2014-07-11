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

double higgs_decay_example(NE6SSM<Two_scale>& ne6ssm, bool speak) {
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
   Eigen::Matrix<double,5,5> Uhiggs_pole;
   int A = 2;  //index for lightest physical pseudo scalar
   //(0 and 1 are goldstones)
   mhiggs = ne6ssm.get_Mhh();
   mAhiggs = ne6ssm.get_MAh();
   Uhiggs = ne6ssm.get_ZH();
   MHmatrix = ne6ssm.get_mass_matrix_hh();
   //   if(fabs(mAhiggs(0) -  90.882)/ mAhiggs(0)  < 1e-3) A = 1;
   mhiggs_pole = ne6ssm.get_physical().Mhh;
   mAhiggs_pole = ne6ssm.get_physical().MAh;
   Uhiggs_pole =  ne6ssm.get_physical().ZH;
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
 GhphysA1A1 = Uhiggs_pole * GHA1A1;
 if(speak){
 std::cout << "GHA1A1 = " << GHA1A1 << std::endl;
 std::cout << "GhA1A1 = " << (Uhiggs_pole * GHA1A1).transpose() << std::endl;
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

    std::cout << "M_1 = M_1^\\prime (GeV) \t & \t";
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
return BRhA1A1;
 
}

bool scanKappaTKappa(int argc, const char *argv[]){
  
  
  double minKapPr=0.01; double maxKapPr=0.1; int stepsKap =300;
  double stepKapPr = (maxKapPr - minKapPr) / (stepsKap-1);
  double minTK=1; double maxTK=50; int stepsTK =300;
  double stepTK = (maxTK - minTK) / (stepsTK-1);
  for(int i=0; i<stepsKap; i++){
     for(int j=0; j<stepsTK; j++){
         GetSpec gs;
         NE6SSM<Two_scale> ne6ssm;
         // NE6SSM<Two_scale> ne6ssm;
         NE6SSM_slha_io slha;
         NE6SSM_input_parameters input_pars;
         int read = gs.input_slha(argc, argv,input_pars);
         if(read != 0) {
            std::cerr << "Fail reading SLHA file." <<std::endl;
            return false;
         }
         input_pars.KappaPrInput = minKapPr + i * stepKapPr;
         input_pars.TKappaPrInput = minTK + j * stepTK;
         int exit_code = gs.getSpectrum(ne6ssm, input_pars);
         if(exit_code != 0) continue;
         double BRhA1A1 = higgs_decay_example(ne6ssm,false);
         const NE6SSM_physical& pole_masses = ne6ssm.get_physical();
         
         const double mHu2 =  ne6ssm.get_mHu2();
         const double yt =  ne6ssm.get_Yu(2,2);
         const double at =  ne6ssm.get_TYu(2,2);
         const double mQ3sq =  ne6ssm.get_mq2(2,2);
         const double mU3sq =  ne6ssm.get_mu2(2,2);
         //this must be positive to evade colour or cahrge breaking minima
         const double ccbfac =  Sqr(yt)*(mQ3sq + mU3sq + mHu2)- Sqr(at);
         
         cout << "  "
              << std::setw(12) << std::left <<ne6ssm.get_KappaPr() << ' '
              << std::setw(12) << std::left << ne6ssm.get_TKappaPr()  << ' '
              << std::setw(12) << std::left << ccbfac << ' '
              << std::setw(12) << std::left << mHu2 << ' '
              << std::setw(12) << std::left << pole_masses.MAh(2) << ' '
              << std::setw(12) << std::left <<  pole_masses.Mhh(0)<< ' '
              << std::setw(12) << std::left << BRhA1A1 << std::endl;
         
         }
  }
  return true;
}


bool scanMQ3MU3At(int argc, const char *argv[]){
   
   
   double minMQ3sq=1.0e+05; double maxMQ3sq=1.0e+08; int stepsMQ3sq =30;
   double stepMQ3sq = log(maxMQ3sq - minMQ3sq) / (stepsMQ3sq-1);
   double minMU3sq=1.0e+05; double maxMU3sq=1.0e+08; int stepsMU3sq =10;
   double stepMU3sq = log(maxMU3sq - minMU3sq) / (stepsMU3sq-1);
   double minTYu=1000; double maxTYu= 3000; int stepsTYu =200;
   double stepTYu = (maxTYu - minTYu) / (stepsTYu-1);
   ofstream mHu2points;
   mHu2points.open("mHu2points.dat");
   for(int i=0; i<stepsMQ3sq; i++){
      for(int j=0; j<stepsTYu; j++){
         for(int k=0; k<stepsMU3sq; k++){
            GetSpec gs;
            NE6SSM<Two_scale> ne6ssm;
            // NE6SSM<Two_scale> ne6ssm;
            NE6SSM_slha_io slha;
            NE6SSM_input_parameters input_pars;
            int read = gs.input_slha(argc, argv,input_pars);
            if(read != 0) {
               std::cerr << "Fail reading SLHA file." <<std::endl;
               return false;
            }
            gs.findSpectrum2(argc,argv,ne6ssm, slha);
            double mHu2 =  ne6ssm.get_mHu2();
            double yt =  ne6ssm.get_Yu(2,2);


            input_pars.mq2Input(2,2) = exp(log(minMQ3sq) + i * stepMQ3sq);
            input_pars.TYuInput(2,2) = minTYu + j * stepTYu;
            input_pars.mu2Input(2,2) = exp(log(minMU3sq) + i * stepMU3sq);
            const double at =  input_pars.TYuInput(2,2);
            const double mQ3sq = input_pars.mq2Input(2,2) ;
            const double mU3sq =input_pars.mu2Input(2,2)  ;
            //this must be positive to evade colour or cahrge breaking minima
            //use 1e+06 for mHu2 to be conservative here.
            double ccbfac =  Sqr(yt)*(mQ3sq + mU3sq+mHu2)- Sqr(at);
            if(ccbfac < 0) continue;
            int exit_code = gs.getSpectrum(ne6ssm, input_pars);
            if(exit_code != 0) continue;
            double BRhA1A1 = higgs_decay_example(ne6ssm,false);
            const NE6SSM_physical& pole_masses = ne6ssm.get_physical();
            //update with correct yt and mHu2
            mHu2 =  ne6ssm.get_mHu2();
            yt =  ne6ssm.get_Yu(2,2);
            ccbfac =  Sqr(yt)*(mQ3sq + mU3sq+mHu2)- Sqr(at);
            
            if(mHu2 < 0){    
              mHu2points << "  "
                         << std::setw(12) << std::left <<mQ3sq << ' ' 
                         << std::setw(12) << std::left <<mU3sq << ' '
                         << std::setw(12) << std::left << at  << ' '
                         << std::setw(12) << std::left << ccbfac << ' '
                         << std::setw(12) << std::left << mHu2 << ' '
                         << std::setw(12) << std::left << pole_masses.MAh(2) << ' '
                         << std::setw(12) << std::left <<  pole_masses.Mhh(0)<< ' '
                         << std::setw(12) << std::left << BRhA1A1 << std::endl;
           }
           else{
              cout << "  "   
                   << std::setw(12) << std::left <<mQ3sq << ' ' 
                   << std::setw(12) << std::left <<mU3sq << ' '
                   << std::setw(12) << std::left << at  << ' '
                   << std::setw(12) << std::left << ccbfac << ' '
                   << std::setw(12) << std::left << mHu2 << ' '
                   << std::setw(12) << std::left << pole_masses.MAh(2) << ' '
                   << std::setw(12) << std::left <<  pole_masses.Mhh(0)<< ' '
                   << std::setw(12) << std::left << BRhA1A1 << std::endl;
           }
      
         }
      }
   }
   mHu2points.close();
   return true;
}


int main(int argc, const char *argv[])
{
   bool speak = true;
   if(speak){
   INFO("=============================");
   INFO("running higgs_decay_example()");
   INFO("=============================");
   }
   // GetSpec gs;
   // NE6SSM<Two_scale> ne6ssm;
   // // NE6SSM<Two_scale> ne6ssm;
   // NE6SSM_slha_io slha;

   // gs.findSpectrum2(argc,argv,ne6ssm, slha);
   
   // higgs_decay_example(ne6ssm,speak);
   
   // if(speak) std::cout << "End of calculation."  << std::endl;
   
   scanMQ3MU3At(argc, argv);

   return 0;
}
