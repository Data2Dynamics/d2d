% This is an example demonstrating how to estimate chain lengths of linear
% chains from data using the identifiability criterion. In the first step,
% data is simulated from a model with known chain length. Then, test models
% with different chain lengths are fitted to this data and the profile
% likelihood is calculated. Using the identifiability criterion, it is then
% possible to find the optimal chain length.

%% SIMULATING DATA
simuCL = 3; %set true chain length
tpoints = [1:3:151]; %set timepoints for simulated data

arInit;
arLoadModel(['LinearChain_' num2str(simuCL)]);
arLoadData('empty')
arCompileAll;

ar.p = [-1 0.5 -1.8 -1.3 -1.3]; %set simulation parameter values
arSimuData(1, 1, tpoints);

filename = sprintf('%s_cl%02d', 'simuData_LinearChain', simuCL);
writeSimuDataToFile(filename, ar, tpoints, simuCL, 1)
fprintf('\n\n*** Simulation completed ***\n\n')


%% RECOVER CHAIN LENGTH FROM SIMULATED DATA
clear ar;

clrange = 1:5; % set range of chain lengths to analyze
for icl = clrange
    arInit;
    arLoadModel(['LinearChain_' num2str(icl)]);
    arLoadData(filename) % load simulated dataset
    arCompileAll
    
    arFitLHS(20); % fit the data with multistart optimization
    
    arPlot 
    arPlotFits  % check if fit is good
    pause(3)
    
    arSave(['LinearChain_cl_' num2str(icl)])
    arPLEInit;
    ple([1 2]) % calculate profiles for k_delay and k_skip
    savefig(['linearChain_profiles_cl' num2str(icl)])
    close all
end

% To find the optimal chain lengths, one has to inspect the profile
% likelihood of k_skip and k_delay (profiles_1.fig, etc.). 
%
% Identifiability criterion: n is the optimal chain length if n+1 is the lowest chain
% length for which k_skip is non-identifiable. If k_delay is non-identifiable for any 
% chain length <=  n, the optimal chain length is 0 instead. If k_skip does not become
% non-identifiable for any chain length, no statement is possible.

function writeSimuDataToFile(filename, ar, tpoints, simuCL, i)
% write .csv file
data = array2table([tpoints' ar.model.data.yExp]);
data.Properties.VariableNames = ['time' ar.model.y];
writetable(data, ['Data/' filename '.csv'],'Encoding','UTF-8');

% write .def file
fileID = fopen(['Data/' filename '.def'], 'w');
fprintf(fileID, ['DESCRIPTION \n\ndata generated from model \n%s'...
    '\n\nwith parameter values\n'], ar.model.description{1});
fprintf(fileID, '%6.6f\n', ar.p);
fprintf(fileID, ['\n\nPREDICTOR \nt T min time 0 100 \n' ...
    '\nINPUTS \n \nOBSERVABLES \n \nERRORS \n \nCONDITIONS \n']);
fclose(fileID);
end
