#include "GetSpec.hpp"

using namespace flexiblesusy;
  
int GetSpec::findSpectrum(int argc, const char* argv[], 
                            NE6SSM<Two_scale> & model, 
                            NE6SSM_slha_io  & slha_io) {


   Command_line_options options(argc, argv);
   if (options.must_print_model_info())
      NE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string rgflow_file(options.get_rgflow_file());
   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   const std::string spectrum_file(options.get_spectrum_file());
   //NE6SSM_slha_io slha_io;
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
   cout << "# "
        << std::setw(12) << std::left << "KapPr" << ' '
        << std::setw(12) << std::left << "TKapPr" << ' '
        << std::setw(12) << std::left << "Mhh(1)" << ' '
        << std::setw(12) << std::left << "MAh(2)" << ' '
        << std::setw(12) << std::left << "BR(h1-> A1 A1)"
        << '\n';
   for(int i = 0; i<= 1000; i++ ){
      for(int j = 1; j<= 100; j++ ){
         input.KappaPrInput = 5e-5 + i*1e-4;
         input.TKappaPrInput = 1 + j;
   
  

   oneset.toMz(); // run SM fermion masses to MZ

   //   NE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
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

   model = spectrum_generator.get_model();
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

   const int exit_code = spectrum_generator.get_exit_code();

   const NE6SSM_physical& pole_masses = model.get_physical();
   const double higgs = pole_masses.Mhh(0);
   const double mA = pole_masses.MAh(2);
   const double mHu2 = model.get_mHu2();
   const double mHd2 = model.get_mHd2();
   const double yt = model.get_Yu(2,2);
   const double at = model.get_TYu(2,2);
   const double mQ3sq = model.get_mq2(2,2);
   const double mU3sq = model.get_mu2(2,2);
    
   // //this must be positive to evade colour or cahrge breaking minima
   const double ccbfac =  Sqr(yt)*(mQ3sq + mU3sq + mHu2)- Sqr(at);
   cout << "  "
        << std::setw(12) << std::left <<model.get_KappaPr() << ' '
        << std::setw(12) << std::left << model.get_TKappaPr()  << ' '
        << std::setw(12) << std::left << ccbfac << ' '
        << std::setw(12) << std::left << mHu2 << ' '
        << std::setw(12) << std::left << mA << ' '
        << std::setw(12) << std::left << higgs;

      }
   }
   return 0;

}
