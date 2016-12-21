% Model of Insulin binding of primary hepatocytes.
% The model is published as part of a multi-scale model in the following
% paper:
% Lars Ole Schwen, Arne Schenk, Clemens Kreutz · María Matilde Bartolomé
% Rodríguez, Lars Kuepfer, Tobias Preusser, Representative Sinusoids For
% Hepatic Four-Scale Pharmacokinetics Simulation, Journal of Mathematical
% Biology, 2014.

%% 
arInit;
arLoadModel('Kreutz_IR_binding');
arLoadData('FacsData_unlog10', 1, 'xls');         
arLoadData('Elisa_relative', 1, 'xls');         

arCompileAll;

%%
global ar

arSetPars('IR_obs_std', []  ,[],[],[], -1.25);
arSetPars('offset_nExpID1', [],[],[],[], 5);
arSetPars('offset_nExpID2', [],[],[],[], 5);
arSetPars('offset_nExpID3', [],[],[],[], 5);
arSetPars('offset_nExpID4', [],[],[],[], 5);
arSetPars('km_nExpID1', [],[],[],3, 8);
arSetPars('km_nExpID2', [],[],[],3, 8);
arSetPars('km_nExpID3', [],[],[],3, 8);
arSetPars('km_nExpID4', [],[],[],3, 8);
arSetPars('scaleElisa_nExpID1', [],[],0,.1,1);
arSetPars('scaleElisa_nExpID2', [],[],0,.1,1);
arSetPars('scaleElisa_nExpID3', [],[],0,.1,1);
arSetPars('scaleElisa_nExpID4', [],[],0,.1,1);
arSetPars('fragments', 1, [], 0, 0, 1)

load Data/prior.mat
fn = fieldnames(prior.mean);
for f=1:length(fn)
    arSetPars(fn{f}, [], [], [], [], [], 1, prior.mean.(fn{f}), prior.std.(fn{f})); % Gaussian prior, sigma(prior) = 3 orders of magnitude
end


