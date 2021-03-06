set terminal x11

set title "NE6SSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='NE6SSM_rgflow.dat'

plot for [i=2:167+1] filename using 1:(column(i)) title columnhead(i)
