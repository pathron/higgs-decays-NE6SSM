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

// File generated at Wed 9 Jul 2014 14:07:20

#include "NE6SSM_input_parameters.hpp"
#include "NE6SSM_slha_io.hpp"
#include "NE6SSM_spectrum_generator.hpp"

#include "error.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      NE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   NE6SSM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;
   QedQcd oneset;
   NE6SSM_input_parameters input;

   if (slha_input_file.empty()) {
      ERROR("No SLHA input file given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   try {
      slha_io.read_from_file(slha_input_file);
      slha_io.fill(oneset);
      slha_io.fill(input);
      slha_io.fill(spectrum_generator_settings);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }
   std::cout << "After reading in from SLHA file we now do scan." << std::endl;
   // const std::vector<double> range (float_range(0., 100., 10));

   cout << "# "
        << std::setw(12) << std::left << "mU3sq" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "MAh(2)"
        << '\n';
   double StartTB = 1;  //double endTB = 50;
   for(int i = 0; i<= 10; i++ ){
      for(int j = 0; j<= 20; j++ ){
         for(int k = 0; k<= 10; k++ ){
   // for (std::vector<double>::const_iterator it = range.begin(),
   //         end = range.end(); it != end; ++it) {
      
      //input.TanBeta = StartTB + i;
      // input.TYuInput(2,2) = -1.73e+03 - i*1e+01;
      // std::cout << "TYu YAMA is "   <<input.TYuInput << std::endl;
            double at =-1.0e+03 - i*2.0e+02; 
            double mU3sq =5e+05 + j*5e+05; 
            double mQ3sq = 5e+05 + k*1e+06; 
            input.mu2Input(2,2)=mU3sq;
            input.mq2Input(2,2)=mQ3sq;
            input.TYuInput(2,2)=at;
            
      oneset.toMz(); // run SM fermion masses to MZ

   NE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   spectrum_generator.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_input_scale(
      slha_io.get_input_scale());
   spectrum_generator.set_parameter_output_scale(
      slha_io.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order));
   spectrum_generator.set_beta_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::beta_loop_order));
   spectrum_generator.set_threshold_corrections(
      spectrum_generator_settings.get(Spectrum_generator_settings::threshold_corrections));

   spectrum_generator.run(oneset, input);

   const NE6SSM<algorithm_type>& model
      = spectrum_generator.get_model();
   const Problems<NE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();

   // output
   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(oneset);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!problems.have_serious_problem())
      slha_io.set_spectrum(model);

   // if (slha_output_file.empty()) {
   //    slha_io.write_to_stream(std::cout);
   // } else {
   //    slha_io.write_to_file(slha_output_file);
   // }

   // if (!spectrum_file.empty())
   //    spectrum_generator.write_spectrum(spectrum_file);

   // if (!rgflow_file.empty())
   //    spectrum_generator.write_running_couplings(rgflow_file);
   const NE6SSM_physical& pole_masses = model.get_physical();
   const double higgs = pole_masses.Mhh(0);
   const double mA = pole_masses.MAh(2);
   const double mHu2 = model.get_mHu2();
   const double mHd2 = model.get_mHd2();
   const double yt = model.get_Yu(2,2);
   //this must be positive to evade colour or cahrge breaking minima
   const double ccbfac = Sqr(at) - Sqr(yt)*(mQ3sq + mU3sq + mHu2);
   const int exit_code = spectrum_generator.get_exit_code();
    cout << "  "
         << std::setw(12) << std::left << mU3sq << ' '
        << std::setw(12) << std::left << mQ3sq << ' '
        << std::setw(12) << std::left << at << ' '
        << std::setw(12) << std::left << ccbfac << ' '
         << std::setw(12) << std::left << mA << ' '
         << std::setw(12) << std::left << higgs << endl;
       
         }
      }   
   }
   //   return exit_codo;
   return 0;
}
