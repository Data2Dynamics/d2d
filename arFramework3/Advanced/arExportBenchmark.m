%Calculate chi2 value without Bessel correction
ar.config.useFitErrorCorrection = false;
if(contains(pwd,'Chen'))
    arSimu(false,false,true);
else
    arFit
    arCalcMerit
end
arWrite_Benchmark
for i = 1:length(ar.model)
    
    for j = 1:length(ar.model(i).data)
        arExportSBML_benchmark(i,j,1);       
    end
end
arWrite_CondXLS
