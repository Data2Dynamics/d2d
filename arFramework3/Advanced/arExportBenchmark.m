%Calculate chi2 value without Bessel correction
ar.config.useFitErrorCorrection = false;
arCalcMerit  
for i = 1:length(ar.model)
    arWrite_Benchmark(i)
    for j = 1:length(ar.model(i).data)
        arExportSBML_benchmark(i,ar.model(i).data(j).cLink,1,j);       
    end
end
arWrite_CondXLS
