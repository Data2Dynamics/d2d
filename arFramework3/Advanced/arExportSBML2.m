% arExportSBML2
% 
% Exports current model to SBML
% arExportSBML_benchmark

function arExportSBML2
    global ar

    %Calculate chi2 value without Bessel correction
    ar.config.useFitErrorCorrection = false;
    if(contains(pwd,'Chen'))
        arCalcMerit
    else
        arFit
        arCalcMerit
    end

    %This file writes out the General info and Model-specific XLS files
    %arWrite_Benchmark

    for i = 1:length(ar.model)   
        for j = 1:length(ar.model(i).data)
    %           Export SBML files
            arExportSBML_benchmark(i,j,1);       
        end
    end