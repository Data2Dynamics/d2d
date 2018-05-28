%%
% addpath('E:\clemens\systemBiologie_cvs\Programmieren\Matlab\matlab-pathway\Library\libSBML-5.12.0-matlab')
% arImportSBML('BIOMD0000000048',100)
arImportSBML('BIOMD0000000048_d2d','tend',100)

%%
Setup

%% Comparison with Biobase:
ar.config.atol=1e-14;
ar.config.rtol=1e-14;

arQplot('x')
arPlot

arCompareWithBiobaseSimulation('SIMU1447940350922.dat');
