for id = 1:length(ar.model(1).data)
    XLS_A = {'time'};
    for iy = 1:length(ar.model(1).data(id).yNames)
        XLS_A = [XLS_A, ar.model(1).data(id).yNames(iy), [ar.model(1).data(id).yNames{iy} '_std']];
    end
    C = ar.model(1).data(id).yExp(:,[1;1]*(1:size(ar.model(1).data(id).yExp,2)));
    C(:,2:2:end) = ar.model(1).data(id).ystdExpSimu;
    T = array2table([ar.model(1).data(id).tExp, C], 'VariableNames',XLS_A);
    out_name = ['./Benchmark_paper/' ar.model(1).name '_condition' num2str(id) '.csv'];
    
    if(~exist('./Benchmark_paper', 'dir'))
        mkdir('./Benchmark_paper')
    end
    
    writetable(T,out_name);
end