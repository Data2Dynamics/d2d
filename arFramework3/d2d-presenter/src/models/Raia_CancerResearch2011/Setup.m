
% Load models & data
arInit;
arLoadModel('il13_jak2_stat5');
arLoadData('MedB1_real_data');
arCompileAll;

% load best fit parameter values
arLoadPars('BestFit');

%arPlot;
arPrint;
