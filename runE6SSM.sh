#!/bin/bash

# Shell script to run the E6SSM in Mathematica, because it
# takes so long.

math<<EOF > sarah_NE6SSM_out.txt
Needs["SARAH`"];
workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = FileNameJoin[{workingDirectory, "Output"}];
SARAH`SARAH[InputDirectories] = {
    ToFileName[{$sarahDir, "Models"}],
    FileNameJoin[{workingDirectory, "sarah"}]
}
Start["NE6SSM"];
ModelOutput[EWSB];
CalcLoopCorrections[EWSB];
CalcRGEs[];
MakeTeX[];
Quit[0];

EOF
