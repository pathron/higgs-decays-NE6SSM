#!/bin/bash 
pdflatex NevzorovSSM-EWSB.tex 
cd Diagrams 
find . -name "*.mp" -exec mpost {} \; 
cd .. 
pdflatex NevzorovSSM-EWSB.tex 
echo "" 
echo "PDF for Model finished" 
echo "Thanks for using SARAH" 
echo "" 
