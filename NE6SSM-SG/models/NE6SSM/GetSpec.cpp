#include "GetSpec.hpp"

using namespace flexiblesusy;

void GetSpec::Output_slha(std::ostream& ostr) {
  
   const Problems<NE6SSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();
   
   const std::string slha_output_file(cmd_line_options.get_slha_output_file());
   const std::string spectrum_file(cmd_line_options.get_spectrum_file());
   const std::string rgflow_file(cmd_line_options.get_rgflow_file());
   // output
   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(oneset);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!problems.have_serious_problem())
      slha_io.set_spectrum(spectrum_generator.get_model());

   if (slha_output_file.empty()) {
      slha_io.write_to_stream(std::cout);
   } else {
      slha_io.write_to_file(slha_output_file);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   return;
}


void GetSpec::Output_slha(std::ostream& ostr, NE6SSM_input_parameters input_pars) {
   //set internal input to match user adjusted one
   input = input_pars;
   //call routine to output to slha file using 
   //private member input 
   Output_slha(std::cout);
   return;

}

int GetSpec::input_slha(int argc, const char* argv[]) {
   ///make sure we start off with empty oneset
   QedQcd newoneset;
   oneset=newoneset;
   //read in slha
   Command_line_options options(argc, argv);
   cmd_line_options = options;
   if (options.must_print_model_info())
      NE6SSM_info::print(std::cout);
   if (options.must_exit())
      return options.status();

   const std::string slha_input_file(options.get_slha_input_file());
   

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
   return 0;
}
int GetSpec::input_slha(int argc, const char* argv[],  NE6SSM_input_parameters & input_pars) {
   //fill priate members from slha file
   int exit_code = input_slha(argc,argv);
   //set user passed input parameters to private member that was
   //filled by slha
   input_pars = input;
   
   return exit_code;
}


int GetSpec::getSpectrum(NE6SSM<Two_scale> & model) {
  
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
   //set model to pass to user
   model = spectrum_generator.get_model();

   const int exit_code = spectrum_generator.get_exit_code();
  
   return exit_code;

}

int GetSpec::getSpectrum(NE6SSM<Two_scale> & model, 
                         NE6SSM_input_parameters  input_pars) {
   input = input_pars;
   getSpectrum(model);

}

int GetSpec::findSpectrum2(int argc, const char* argv[], 
                            NE6SSM<Two_scale> & model, 
                            NE6SSM_slha_io  & slha_io) { 
   int read_ok = input_slha(argc, argv);
   if(read_ok!=0) return read_ok;
   int exit_code = getSpectrum(model);
   Output_slha(std::cout);
   
   return exit_code;

}
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

   if (slha_output_file.empty()) {
      slha_io.write_to_stream(std::cout);
   } else {
      slha_io.write_to_file(slha_output_file);
   }

   if (!spectrum_file.empty())
      spectrum_generator.write_spectrum(spectrum_file);

   if (!rgflow_file.empty())
      spectrum_generator.write_running_couplings(rgflow_file);

   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;

}
