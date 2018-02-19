function arWrite_CondXLS(imodel)
global ar
if(~exist('imodel','var') || isempty(imodel))
    imodel = 1:length(ar.model);
end

for im = imodel
%     %old code to include conditions in data file
%     tmp_par = {};
%     for jc= 1:length(ar.model(im).condition)
%             tmp_id = ar.model(im).condition(jc).dLink(1);
% 
%             for jd=1:length(ar.model(im).data(tmp_id).condition)
%                 if(sum(strcmp(ar.model(im).data(tmp_id).condition(jd).parameter,tmp_par))==0)
%                     tmp_par = [tmp_par ar.model(im).data(tmp_id).condition(jd).parameter];
%                 end
%             end
%     end

    for id = 1:length(ar.model(im).data)
        XLS_Data = {'time'};%, tmp_par{:}};
        %tmp_C = NaN(length(ar.model(im).data(id).tExp),length(tmp_par));

        for iy = 1:length(ar.model(im).data(id).yNames)
            XLS_Data = [XLS_Data, ar.model(im).data(id).yNames(iy), [ar.model(im).data(id).yNames{iy} '_std']];           
        end
        XLS_Sim = [{'time'}, ar.model(im).data(id).yNames];
        Exp_Data = [ar.model(im).data(id).tExp ar.model(im).data(id).yExp(:,[1;1]*(1:size(ar.model(im).data(id).yExp,2)))];
        Exp_Data(:,3:2:end) = ar.model(im).data(id).ystdExpSimu;
        Sim_Data = [ar.model(im).data(id).tExp ar.model(im).data(id).yExpSimu];
        %old code to include conditions in data file
%         C = [tmp_C C];
%         [tmp,ia] = setdiff(ar.model(im).data(id).pold,ar.model(im).data(id).fp);
%         ia(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
%         tmp(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
%         for i = 1:length(tmp)
%            cond_id = find(ismember(tmp_par,tmp{i}));
%            if(~isempty(cond_id) && length(cond_id)==1)
%             C(:,cond_id) = repmat(str2double(strrep(strrep(ar.model(im).data(id).fp{ia(i)},'(',''),')','')),length(ar.model(im).data(id).tExp),1);
%            end
%         end
%         
%         nr_of_condition = find(ar.model(im).condition(ar.model(im).data(id).cLink).dLink==id);

        out_name = ['./Benchmark_paper/Data/' ar.model(im).name '_m' num2str(im) '_data' num2str(id) '.xlsx'];

        if(~exist('./Benchmark_paper', 'dir'))
            mkdir('./Benchmark_paper')
        end
        if(~exist('./Benchmark_paper/Data', 'dir'))
            mkdir('./Benchmark_paper/Data')
        end
        
        %prepare excel file
        xlwrite(out_name,XLS_Data,'Exp Data');
        xlwrite(out_name,Exp_Data,'Exp Data','A2');
        xlwrite(out_name,XLS_Sim,'Simulation');
        xlwrite(out_name,Sim_Data,'Simulation','A2');
    end
end