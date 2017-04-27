function arWrite_CondXLS(imodel)
global ar
if(~exist('imodel','var') || isempty(imodel))
    imodel = 1:length(ar.model);
end

for im = imodel
    tmp_par = {};
    for jc= 1:length(ar.model(im).condition)
            tmp_id = ar.model(im).condition(jc).dLink(1);

            for jd=1:length(ar.model(im).data(tmp_id).condition)
                if(sum(strcmp(ar.model(im).data(tmp_id).condition(jd).parameter,tmp_par))==0)
                    tmp_par = [tmp_par ar.model(im).data(tmp_id).condition(jd).parameter];
                end
            end
    end

    for id = 1:length(ar.model(im).data)
        XLS_A = {'time', tmp_par{:}};
        tmp_C = NaN(length(ar.model(im).data(id).tExp),length(tmp_par));

        for iy = 1:length(ar.model(im).data(id).yNames)
            XLS_A = [XLS_A, ar.model(im).data(id).yNames(iy), [ar.model(im).data(id).yNames{iy} '_std']];
        end
        C = ar.model(im).data(id).yExp(:,[1;1]*(1:size(ar.model(im).data(id).yExp,2)));
        C(:,2:2:end) = ar.model(im).data(id).ystdExpSimu;
        C = [tmp_C C];

        [tmp,ia] = setdiff(ar.model(im).data(id).pold,ar.model(im).data(id).fp);
        ia(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
        tmp(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
        for i = 1:length(tmp)
           cond_id = find(ismember(tmp_par,tmp{i}));
           if(~isempty(cond_id) && length(cond_id)==1)
            C(:,cond_id) = repmat(str2double(strrep(strrep(ar.model(im).data(id).fp{ia(i)},'(',''),')','')),length(ar.model(im).data(id).tExp),1);
           end
        end

        T = array2table([ar.model(im).data(id).tExp, C], 'VariableNames',XLS_A);
        out_name = ['./Benchmark_paper/' ar.model(im).name '_m' num2str(im) '_condition' num2str(id) '.csv'];

        if(~exist('./Benchmark_paper', 'dir'))
            mkdir('./Benchmark_paper')
        end

        writetable(T,out_name);
    end
end