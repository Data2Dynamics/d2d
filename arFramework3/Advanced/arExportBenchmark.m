%Calculate chi2 value without Bessel correction
ar.config.useFitErrorCorrection = false;
ar.config.fiterrors = 1;
if(contains(pwd,'Chen'))
    arSimu(false,false,true);
else
    arFit
    arCalcMerit
end

%This file writes out the General info and Model-specific XLS files
arWrite_Benchmark

% for i = 1:length(ar.model)   
%     for j = 1:length(ar.model(i).data)
%           Export SBML files
%         arExportSBML_benchmark(i,j,1);       
%     end
% end

%Export data values, uncertainties and model simulations
% arWrite_CondXLS
