arCalcMerit
for i = 1:length(ar.model)
    arWrite_Benchmark(i)
end
arWrite_CondXLS
arExportSBML(1,1,0)