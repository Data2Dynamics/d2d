% File to export human-readable formats of a model
% For every data entry in the struct, the measurements are written out in a
% single file, together with the model value at these time points
% 
% Regarding the data, experimental errors as well as the uncertainties
% estimated via error model are stored
%
% For the corresponding ODE and parameter output, see arWrite_Benchmark

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
            XLS_Data = [XLS_Data, ar.model(im).data(id).yNames(iy), [ar.model(im).data(id).yNames{iy} '_std'], [ar.model(im).data(id).yNames{iy} '_ExpError']];           
        end
        XLS_Sim = [{'time'}, ar.model(im).data(id).yNames];
        Exp_Data = [ar.model(im).data(id).tExp ar.model(im).data(id).yExp(:,[1;1;1]*(1:size(ar.model(im).data(id).yExp,2)))];
        
        %Prefer Exp Errors, only set model errors for empty ones
        Exp_Data(:,3:3:end) = NaN;
        Exp_Data(:,4:3:end) = ar.model(im).data(id).yExpStd;
        
        tmp_stdSimu = ar.model(im).data(id).ystdExpSimu;
        tmp_stdSimu(~isnan(ar.model(im).data(id).yExpStd)) = NaN;
        Exp_Data(:,3:3:end) = tmp_stdSimu;
        Sim_Data = [ar.model(im).data(id).tExp ar.model(im).data(id).yExpSimu];

        out_name = ['./Benchmark_paper/Data/' 'model' num2str(im) '_data' num2str(id) '.xlsx'];

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
        
%         %Prepare for simu without error model, copy estimated uncertainties
        ar.model(im).data(id).yExpStd_cp = ar.model(im).data(id).yExpStd;
        ar.model(im).data(id).yExpStd(isnan(ar.model(im).data(id).yExpStd)) = ar.model(im).data(id).ystdExpSimu(isnan(ar.model(im).data(id).yExpStd));
    end
end 

%Don't do fitting without error model, and separate sheets for pars/simulations
return;

    p_backup = ar.p;
    %Do simulations without error model
    bkp_errorPars = ar.qFit(~cellfun(@isempty,strfind(ar.pLabel,'sd_')));
    bkp_fiterror = ar.config.fiterrors;
    arSetParsPattern('sd_',[],0)
    ar.config.fiterrors=-1;
    
    if(contains(pwd,'Chen'))
    	arSimu(false,false,true);
    else
        arFit
    end
    chi2Val = arGetMerit;
    
    Chi2_string = {'Chi2 without Error model',chi2Val};    
    xlwrite(['./Benchmark_paper/General_info.xlsx'],Chi2_string,'General Info','A8:B8');
    
%     return;
    % Write parameter values to the general info file
    General_string(1,1:6) = {'parameter','value','lower boundary','upper boundary','analysis at log-scale','estimated'};    
    
    for i = find(~ar.qError)
        if ar.qLog10(i)
            General_string{end+1,1} = ['log10(' ar.pLabel{i} ')']; 
        else
            General_string{end+1,1} = ar.pLabel{i}; 
        end
        General_string{end,2} = ar.p(i);
        General_string{end,3} = ar.lb(i);
        General_string{end,4} = ar.ub(i);
        General_string{end,5} = ar.qLog10(i);
        if(ar.qFit(i) == 1)
            General_string{end,6} = 'yes';
        else
            General_string{end,6} = 'fixed';
        end
    end   
    xlwrite(['./Benchmark_paper/General_info.xlsx'],General_string,'Parameters_noErrorModel');
    
for im = imodel    
    for id = 1:length(ar.model(im).data)
        out_name = ['./Benchmark_paper/Data/' 'model' num2str(im) '_data' num2str(id) '.xlsx'];
        XLS_Sim = [{'time'}, ar.model(im).data(id).yNames];
        Sim_Data = [ar.model(im).data(id).tExp ar.model(im).data(id).yExpSimu];

        xlwrite(out_name,XLS_Sim,'Simulation_noErrorModel');
        xlwrite(out_name,Sim_Data,'Simulation_noErrorModel','A2');
        
        ar.model(im).data(id).yExpStd = ar.model(im).data(id).yExpStd_cp;
    end
end
%Reset error model
ar.config.fiterrors=bkp_fiterror;
ar.qFit(~cellfun(@isempty,strfind(ar.pLabel,'sd_'))) = bkp_errorPars;
ar.p = p_backup;