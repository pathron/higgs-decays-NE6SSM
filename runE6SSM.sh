#!/bin/bash

# Shell script to run the E6SSM in Mathematica, because it
# takes so long.

math<<EOF > sarah_E6SSM_out.txt

<<SARAH\`;
Start["E6SSM"];
ModelOutput[EWSB];
CalcLoopCorrections[EWSB];
CalcRGEs[];
MakeTeX[];
Quit[0];

EOF
