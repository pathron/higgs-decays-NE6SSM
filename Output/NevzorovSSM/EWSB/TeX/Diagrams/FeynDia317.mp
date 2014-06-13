% Diagrams/FeynDia317.mp -- do not edit, generated automatically by NevzorovSSM-EWSB.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 75*1.0pt#, 75*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__l1);
vright (__r1, __r2);
vconnect ("plain", __l1, __v1);
vconnect ("plain", __r1, __v1);
vconnect ("wiggly", __r2, __v1);
vlabel ("$\tilde {\chi }^{'0}_{{i}}$", __l1);
vlabel ("$\tilde {\chi }^{'0}_{{j}}$", __r1);
vlabel ("${Z'}_{{\mu }}$", __r2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
