% File to export human-readable formats of a model
% Generates a general-info file with fact sheet, parameters, raw ODEs,
% parameter transformations (including initials) and experimental conditions
%
% Plus for every data entry in the struct, a file with finalized ODEs, with
% all transformations applied, is written out, together with the initial
% values and all observation functions and possible log transform
%
% For corresponding output of data and model simulations, see
% arWrite_CondXLS

function arWrite_Benchmark
global ar
    if(~exist('./Benchmark_paper', 'dir'))
        mkdir('./Benchmark_paper')
    end
    if(~exist('./Benchmark_paper/Model', 'dir'))
        mkdir('./Benchmark_paper/Model')
    end
    
    %% General csv sheets
    
    %General Infos
    [~,chi2Val] = arGetMerit;
    
    model_nameStr = strsplit(pwd,'/');
    
    %Set Model names manually
    if(strcmp(model_nameStr{end},"Beer_MolBiosyst2014"))
        model_name = "Beer_MolBioSystems2014";
    elseif(strcmp(model_nameStr{end},"TGFb_ComplexModel_WithGenes_Reduced"))
        model_name = "Lucarelli_CellSystems2017";
    elseif(strcmp(model_nameStr{end},"Merkle_JAK2STAT5_PCB2016"))
        model_name = "Merkle_PCB2016";
    elseif(strcmp(model_nameStr{end},"Reelin_PONE2017"))
        model_name = "Hass_PONE2017";
    elseif(strcmp(model_nameStr{end},"Schwen_InsulinMouseHepatocytes_PlosOne2014"))
        model_name = "Schwen_PONE2014";
    elseif(strcmp(model_nameStr{end},"Crauste_ImmuneCells_CellSystems2017"))
        model_name = "Crauste_CellSystems2017";    
    elseif(strcmp(model_nameStr{end},"Bruno_Carotines_JExpBio2016"))
        model_name = "Bruno_JExpBio2016";
    else
        model_name = model_nameStr{end};
    end
    
    General_string = {'Model:'; ''; 'FACTS'};
    General_string{1,2} = model_name;
    sum_cond = 0;
    for im = 1:length(ar.model)
       sum_cond = sum_cond + length(ar.model(im).condition); 
    end
    tmp = sprintf('The model contains %i data points, %i free parameters and %i experimental conditions',ar.ndata,sum(ar.qFit==1),sum_cond);
    General_string{end+1,1} = tmp;
    General_string{end+1,1} = 'Errors are assumend as additive Gaussian errors';
    if ar.config.fiterrors == -1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)==1)==0)
        General_string(end+2,1:2) = {'Chi2 value of the model is ', chi2Val.chi2_all};
        General_string(end+3,1:2)={'Compartments','size'};
    else
        General_string(end+2,1:2) = {'Likelihood value of the model is ', chi2Val.loglik_all};
        General_string(end+1,1:2) = {'Chi2 value of the model is ', chi2Val.chi2_res + chi2Val.chi2_prior + chi2Val.chi2_constr};
        General_string(end+2,1:2)={'Compartments','size'};
    end
    
    for im = 1:length(ar.model)
        for i = 1:length(ar.model(im).c)
            General_string = [General_string; {ar.model(im).c{i}, str2num(ar.model(im).pc{i})}];
        end
    end
    
    General_string{end+2,1} = 'Special features of the model';
    General_string{end+1,1} = 'Splines';
    splinestr = 'no';
    for im = 1:length(ar.model)        
        if(contains(ar.model(im).fu,'spline'))
            splinestr = 'yes';
        end
    end
    General_string{end,2} = splinestr;
    %General_string{end+1,1} = 'Estimated error parameters';
    est_errors(1:length(ar.model(im).data)) = false;
    log_fit = 'no';
    for im = 1:length(ar.model)
        for id = 1:length(ar.model(im).data)
            %if(sum(sum(isnan(ar.model(im).data(id).yExpStd)))>sum(sum(isnan(ar.model(im).data(id).yExp))) && sum(sum(isnan(ar.model(im).data(id).yExp)))<(length(ar.model(im).data(id).tExp)*size(ar.model(im).data(id).yExp,2)))
            if(sum(isnan(ar.model(im).data(id).yExpStd) & ~isnan(ar.model(im).data(id).yExp))>0)
                est_errors(id) = true;    
            end
            if(sum(ar.model(im).data(id).logfitting)>0)
               log_fit = 'yes'; 
            end
        end
    end
    if(all(est_errors==false))
        est_error_string = 'none';
    elseif(any(est_errors == false))
        est_error_string = 'mix of fixed and estimated';
    elseif(all(est_errors==true))
        est_error_string = 'all';
    end
    %General_string{end,2} = est_error_string;
    
    if(strcmp(est_errors,'yes'))
       %General_string{end,3} = 'Check sheet Parameters_noErrorModel for values with fixed errors'; 
    end
    
    General_string{end+1,1} = 'log-scale of observations';
    General_string{end,2} = log_fit;
    
    General_string{end+1,1} = 'Events';
    General_string{end,2} = ar.config.useEvents;
    
    General_string{end+1,1} = 'Number of parameter dependent inital conditions';
    General_string{end,2} = sum(ar.qInitial & ar.qFit);   
        
    
    
    xlwrite(['./Benchmark_paper/General_info.xlsx'],General_string,'General Info');
    
     % Parameter values
    General_string = {};
    General_string(1,1:6) = {'parameter','value','lower boundary','upper boundary','analysis at log-scale','estimated'};
    
    for i = 1:length(ar.p)
        if(ar.qError(i) == 1 && (ar.config.fiterrors == -1 || (ar.p(i) == -1 && ar.qFit(i)~=1) || strcmp(est_error_string,'none')))
            continue;
        end
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
    
    xlwrite(['./Benchmark_paper/General_info.xlsx'],General_string,'Parameters');
    
    % Conditions
    Conditions_string = {'Conditions'};    
    
    
    for im = 1:length(ar.model)
        for jc = 1:length(ar.model(im).condition)
            if(jc==1)
                pold_tmp = ar.model(im).condition(jc).pold;
                continue;
            end
            not_inList = ~ismember(ar.model(im).condition(jc).pold,pold_tmp);
            if(~isempty(not_inList) && sum(not_inList)>0)
                pold_tmp(end+1:end+sum(not_inList)) = ar.model(im).condition(jc).pold(not_inList);
            end
        end
        
        %loop over data to get condition values
        for jd= 1:length(ar.model(im).data)
            cond_tmp = cell(2,length(ar.model(im).data(jd).condition));
            for jc = 1:length(ar.model(im).data(jd).condition)
               cond_tmp{1,jc} =  ar.model(im).data(jd).condition(jc).parameter;
               cond_tmp{2,jc} =  ar.model(im).data(jd).condition(jc).value;
            end
            cond_notIn = ~ismember(cond_tmp(1,:),pold_tmp);
            if(~isempty(cond_notIn) && sum(cond_notIn)>0)
                pold_tmp(end+1:end+sum(cond_notIn)) = cond_tmp(1,cond_notIn);
            end
        end
    
        par_trafo = cell(length(pold_tmp),length(ar.model(im).data));
        num_cols = 6;
        %Collect parameter trafos
        Conditions_string(end+2,1:num_cols) = {['Model file ' num2str(im)],'','exp condition','nTimePoints','nDataPoints','chi2 value'}; 
        for jd= 1:length(ar.model(im).data)
           jc = ar.model(im).data(jd).cLink;          
           Conditions_string{end+1,2} = ['Data file ' num2str(jd)];
           Conditions_string{end,3} = jc;
           Conditions_string{end,4} = length(ar.model(im).data(jd).tExp);
           Conditions_string{end,5} = arnansum(ar.model(im).data(jd).ndata);
           Conditions_string{end,6} = arnansum(ar.model(im).data(jd).chi2);

           %append parameter trafos in each data struct
           par_trafo(ismember(pold_tmp,ar.model(im).condition(jc).pold),jd) = regexprep(regexprep(ar.model(im).condition(jc).fp,'^(',''),')$','');
           
           %Append condition values
           cond_tmp = cell(2,length(ar.model(im).data(jd).condition));
            for jc = 1:length(ar.model(im).data(jd).condition)
               cond_tmp{1,jc} =  ar.model(im).data(jd).condition(jc).parameter;
               cond_tmp{2,jc} =  ar.model(im).data(jd).condition(jc).value;
            end
            [which_cond,which_pold] = ismember(cond_tmp(1,:),pold_tmp);
            par_trafo(which_pold,jd) = cond_tmp(2,which_cond);
            
        end
        
        %get differences in parameter trafos
        num = zeros(1,size(par_trafo,1));
        istrafo = false(1,size(par_trafo,1));
        for i = 1:size(par_trafo,1)
            ids_empty = cellfun(@isempty,par_trafo(i,:));
            if(any(ids_empty))
                for ie = find(ids_empty==1)
                    par_trafo{i,ie} = '';
                end
            end
            num(i) = length(unique(par_trafo(i,:)));
            istrafo(i) = num(i)==1 & ~strcmp(par_trafo{i,1},pold_tmp{i});
        end   
        
      %Comment out if for 1 condition models, all trafos should be printed
%         if(length(ar.model(im).condition)==1)
%             tmp_par = ar.model(im).condition(max_id).pold(num>=1);
%         else
            tmp_par = pold_tmp(num>1);
%         end
        %Go through differing parameter trafos (num>1)
        Conditions_string(end-length(ar.model(im).data),(num_cols+1):length(tmp_par)+num_cols) = tmp_par;
%         if(length(ar.model(im).condition)==1)
%             Conditions_string(end-(length(ar.model(im).data)-1):end,6:length(tmp_par)+5) = par_trafo(num>=1,:)';
%         else
            Conditions_string(end-(length(ar.model(im).data)-1):end,(num_cols+1):length(tmp_par)+num_cols) = par_trafo(num>1,:)';%strrep(strrep(par_trafo(num>1,:)','(',''),')','');
%         end

        if(im == 1)
            Conditions_string_tmp = cell(1,2);
            Conditions_string_tmp(1,1) = {'General transformations'};
            Conditions_string_tmp(end+1,1:2) = {'Parameter','Replacement'};
        end
        if(sum(istrafo)>0)
            Conditions_string_tmp(end+1:size(Conditions_string_tmp,1)+sum(istrafo),1) = pold_tmp(istrafo)';   
            Conditions_string_tmp(size(Conditions_string_tmp,1)-sum(istrafo)+1:size(Conditions_string_tmp,1),2) = par_trafo(istrafo,1);%strrep(strrep(ar.model(im).condition(jc).fp(istrafo),'(',''),')','');        
        end
        for us = 1:length(ar.model(im).u)
            Conditions_string_tmp{end+1,1} = ar.model(im).u{us};
            Conditions_string_tmp{end,2} = ar.model(im).fu{us};
        end
        
    end
    
    %This part writes out the constants (num==1 means same parameter trafo
    %in every condition)
    if(size(Conditions_string_tmp,1)>2)
        Conditions_string(end+3:size(Conditions_string,1)+size(Conditions_string_tmp,1)+2,1:2) = Conditions_string_tmp;
    end
    
    
    xlwrite(['./Benchmark_paper/General_info.xlsx'],Conditions_string,'Experimental conditions');
   
    
    % Write out RAW ODEs, to see parameter trafos and dependencies
    rawODE_string = {'ODE equations'};
    for im = 1:length(ar.model)
        rawODE_string{end+2,1} = ['Model ' num2str(im)];
        tmp_ode = ar.model(im).fx;
        for i = 1:length(ar.model(im).x)
            rawODE_string{end+1,1} = ['d' ar.model(im).x{i} '/dt'];                
            rawODE_string{end,2} = char(sym(tmp_ode{i}));                
        end
    end
    
    xlwrite(['./Benchmark_paper/General_info.xlsx'],rawODE_string,'Raw ODEs');   
    
    
%     return;
    %% Data specific csv files
    
    %ODE equations
    for im = 1:length(ar.model)
        for id = 1:length(ar.model(im).data)
            file_name = ['./Benchmark_paper/Model/' 'model' num2str(im) '_data' num2str(id) '.xlsx'];

            if(~isempty(ar.model(im).x))
                jc = ar.model(im).data(id).cLink;
                %file_name = ['./Benchmark_paper/Model/' 'model' num2str(im) '_condition' num2str(jc) '.xlsx'];
                tmp_ode = ar.model(im).fx;
                for j = 1:length(ar.model(im).u)
                    regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(im).u{j},ar.model(im).u{j},ar.model(im).u{j},ar.model(im).u{j});
                    tmp_ode = regexprep(tmp_ode,regPar, ['(' ar.model(im).fu{j} ')']);    
                end

                ode_string = {'ODE equations'; ''};
                tmp_ode2 = tmp_ode;
                for j = 1:length(ar.model(im).condition(jc).pold)
                    regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(im).condition(jc).pold{j},ar.model(im).condition(jc).pold{j},ar.model(im).condition(jc).pold{j},ar.model(im).condition(jc).pold{j});
                    %tmp_ode2 = regexprep(tmp_ode2,regPar, strrep(strrep(ar.model(im).condition(jc).fp{j},'(',''),')',''));  
                    tmp_ode2 = regexprep(tmp_ode2,regPar, ar.model(im).condition(jc).fp{j});
                end
                for i = 1:length(ar.model(im).x)
                    ode_string{i+2,1} = ['d' ar.model(im).x{i} '/dt'];                
                    ode_string{i+2,2} = char(sym(tmp_ode2{i}));                
                end

    %             % Inputs    
    %             if(~isempty(ar.model(im).fu))
    % 
    %                 ode_string(length(ar.model(im).x)+3:length(ar.model(im).x)+5,1:2) = {'', '';'The following', ' parameters are obtained from input functions';'Parameter', 'Input function'};
    % 
    %                 for i=1:length(ar.model(im).condition(jc).fu)
    %                     ode_string{end+1,1} = ar.model(im).u{i};
    %                     ode_string{end,2} = ar.model(im).condition(jc).fu{i};
    %                 end
    %             end

                xlwrite(file_name,ode_string,'ODEs');

            end  

        %% Observations

            obs_string = {'Observables';''};

            obs_string(2,2:5) = {'scale','observation function','uncertainties','error model'};

            for i = 1:length(ar.model(im).data(id).yNames)            

                if(sum(isnan(ar.model(im).data(id).yExpStd(:,i)))==length(ar.model(im).data(id).tExp) && sum(isnan(ar.model(im).data(id).yExp(:,i)))<length(ar.model(im).data(id).tExp))
                   A_tmp = 'fitted';
                elseif(sum(isnan(ar.model(im).data(id).yExpStd(:,i)))==sum(isnan(ar.model(im).data(id).yExp(:,i))) && sum(isnan(ar.model(im).data(id).yExp(:,i)))<length(ar.model(im).data(id).tExp))
                   A_tmp = 'measured';
                else
                   A_tmp = 'mixed';
                end


                if(sum(ar.model(im).data(id).logfitting(i))==0)
                    tmp_text = 'non-log';
                elseif(sum(ar.model(im).data(id).logfitting(i))==1)
                    tmp_text = 'log';
                end
               obs_string{end+1,1} = ar.model(im).data(id).yNames{i};
               obs_string{end,2} = tmp_text;
               if(strcmp(tmp_text,'log'))
                   obs_string{end,3} = ['log10(' ar.model(im).data(id).fy{i} ')'];
               else
                   obs_string{end,3} = ar.model(im).data(id).fy{i};
               end
               obs_string{end,4} = A_tmp;
               obs_string{end,5} = ar.model(im).data(id).fystd{i};

               %Deprecated, normalization of data to 1
               %obs_string{end,6} = ar.model(im).data(id).normalize(i);
            end


            if(~isempty(ar.model(im).fz))
                found_z = 0;
                obs_string{end+2,1} = '';
                for jz = 1:length(ar.model(im).z)
                    regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(im).z{jz},ar.model(im).z{jz},ar.model(im).z{jz},ar.model(im).z{jz});
                    z_index = regexp(ar.model(im).data(id).fy,regPar);
                    check_empty = cellfun(@isempty,z_index);
                    if(sum(check_empty)<length(check_empty))
                        obs_string{end+1,1} = ar.model(im).z{jz};
                        obs_string{end,2} = ar.model(im).fz{jz};
                        found_z = found_z+1;
                    end
                end  
                if(found_z>0)
                    obs_string{end-found_z,1} = 'With definitions';
                end
            end
            xlwrite(file_name,obs_string,'Observables');

        %% Initials

            init_string = {'Initial values';''};
            tmp = ~cellfun(@isempty,strfind(ar.model(im).condition(jc).pold,'init_'));

            for i = find(tmp)
               init_string{end+1,1} = ar.model(im).condition(jc).pold{i};
               %init_string{end,2} = strrep(strrep(ar.model(im).condition(jc).fp{i},'(',''),')','');
               init_string{end,2} = ar.model(im).condition(jc).fp{i};
            end

            xlwrite(file_name,init_string,'Initials');
        end
    end
    %% Constants
    
%     A = {'Constant values';'Condition 1 is taken as placeholder, check if constants are representative for model!'};
%     [tmp,ia] = setdiff(ar.model(im).condition(1).pold,ar.model(im).condition(1).fp);
%     ia(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
%     tmp(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
%     for i = 1:length(tmp)
%        A{end+1,1} = tmp{i};
%        A{end,2} = strrep(strrep(ar.model(im).condition(1).fp{ia(i)},'(',''),')','');
%     end
%     
%     fid = fopen(['./Benchmark_paper/constants_m' num2str(im) '.csv'],'w');
%     for i = 1:size(A,1)
%         fprintf(fid, '%s; %s \n',A{i,1},A{i,2});
%     end
%     fclose(fid)
end
