tic;
% Load models & data
arInit;
arLoadModel('il13_jak2_stat5');
arLoadData('MedB1_real_data');
arLoadData('MedB1_real_data2');
arCompileAll;

% load best fit parameter values
arLoadPars('BestFit');
toc;
%arPlot;
%arPrint;
