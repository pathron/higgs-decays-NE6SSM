% Diagrams/FeynDia161.mp -- do not edit, generated automatically by NevzorovSSM-EWSB.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 75*1.0pt#, 75*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__l1, __l2);
vright (__r1, __r2);
vconnect ("dashes", __l1, __v1);
vconnect ("dashes", __r1, __v1);
vconnect ("wiggly", __r2, __v1);
vconnect ("wiggly", __l2, __v1);
vlabel ("$h_{{i}}$", __l1);
vlabel ("$h_{{j}}$", __r1);
vlabel ("$Z_{{\mu }}$", __r2);
vlabel ("$Z_{{\nu }}$", __l2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
