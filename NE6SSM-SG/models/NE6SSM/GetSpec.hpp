#ifndef GETSPEC_H
#define GETSPEC_H

#include "NE6SSM_input_parameters.hpp"
#include "NE6SSM_slha_io.hpp"
#include "NE6SSM_spectrum_generator.hpp"
#include "spectrum_generator_settings.hpp"
#include "error.hpp"
#include "lowe.h"
#include "command_line_options.hpp"
#include "NE6SSM_two_scale_model.hpp"
#include "two_scale_model.hpp"

#include <iostream>
#include <cstdlib>

namespace flexiblesusy {

class GetSpec {

private: 
   NE6SSM_spectrum_generator<Two_scale> spectrum_generator;
   Spectrum_generator_settings spectrum_generator_settings;
   NE6SSM_slha_io  slha_io;
   NE6SSM_input_parameters input;
   QedQcd oneset;
   Command_line_options cmd_line_options;
   
public:
   
   ///fills private members slha_io, input and  
   int input_slha(int argc, const char* argv[]);
   ///fills private mebers and sets passed input for public use
   int input_slha(int argc, const char* argv[],NE6SSM_input_parameters & inp);
   ///oututs slha spectra to user specified ostream
   void Output_slha(std::ostream& ostr);
   ///sets input according to user passed object 
   ///oututs slha spectra to user specified ostream 
   void Output_slha(std::ostream& ostr, NE6SSM_input_parameters input_pars);
   
   
   ///finds the spectrum after slha_io and input have been set   
   int getSpectrum(NE6SSM<Two_scale> & model);
   ///finds the spectrum for passed input_pars after spectrum has been read in  
   int getSpectrum(NE6SSM<Two_scale> & model, NE6SSM_input_parameters inp);
   
   ///sets up from slha input file and does full spectrum calculation 
   ///for the user but passes out a NE6SSM<Two_scale> and  NE6SSM_slha_io
   ///object
   int findSpectrum(int argc, const char* argv[], NE6SSM<Two_scale> & model, 
                      NE6SSM_slha_io  & slha);
   
   //as above but done by calling other routines so neater
   int findSpectrum2(int argc, const char* argv[], NE6SSM<Two_scale> & model, 
                     NE6SSM_slha_io  & slha_io);

};


} // namespace flexiblesusy



#endif
