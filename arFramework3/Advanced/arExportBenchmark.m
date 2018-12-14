% arExportBenchmark
% 
% Exports current model for benchmark database, including general info, 
% model XLS, SBMLs, data, uncertainties, model simulations 
% Calls serial: arWrite_Benchmark, arExportSBML_benchmark, arWrite_CondXLS
% (fits the model before exporting)
% 
% See also arWrite_Benchmark, arExportSBML_benchmark, arWrite_CondXLS

function arExportBenchmark

    %Calculate chi2 value without Bessel correction
    ar.config.useFitErrorCorrection = false;
    if(contains(pwd,'Chen'))
        arCalcMerit
    else
        arFit
        arCalcMerit
    end

    %This file writes out the General info and Model-specific XLS files
    arWrite_Benchmark

    for i = 1:length(ar.model)   
        for j = 1:length(ar.model(i).data)
    %           Export SBML files
            arExportSBML_benchmark(i,j,1);       
        end
    end

    % Export data values, uncertainties and model simulations
    arWrite_CondXLS
