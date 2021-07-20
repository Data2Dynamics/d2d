% Run this file once in order to compile the mexbobyqa.F file to the mexbobyqa.mex***
% format, which can be recognized by Matlab. Then you will be able to use the 
% bobyqa.m function.

% Required: Your current Matlab installation must be connected to a Fortran compiler. 
% For further information see: 
% https://mathworks.com/help/matlab/build-fortran-mex-files-1.html

mex mexbobyqa.F
