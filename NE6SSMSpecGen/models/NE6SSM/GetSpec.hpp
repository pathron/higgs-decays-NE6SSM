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

class Two_scale;

namespace flexiblesusy {

class GetSpec {

private: 
   NE6SSM_spectrum_generator<Two_scale> spectrum_generator;

public:
   int findSpectrum(int argc, const char* argv[], NE6SSM<Two_scale> & model, 
                      NE6SSM_slha_io  & slha);
   
};

} // namespace flexiblesusy



#endif
