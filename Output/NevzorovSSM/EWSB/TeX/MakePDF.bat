pdflatex NevzorovSSM-EWSB.tex 
cd Diagrams 
FOR %%I IN (*.mp) DO MPOST "%%I" 
cd .. 
pdflatex NevzorovSSM-EWSB.tex 
