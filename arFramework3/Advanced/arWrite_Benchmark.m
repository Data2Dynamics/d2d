% arWrite_Benchmark
% 
% Export function used in Hass et al, 2018.
% 
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
% 
% See also arCopyBenchmarkModels
% 
% References: 
% Helge Hass, Carolin Loos, Elba Raimundez Alvarez, Jens
% Timmer, Jan Hasenauer, Clemens Kreutz, Benchmark Problems for Dynamic
% Modeling of Intracellular Processes doi: https://doi.org/10.1101/404590  


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
    General_string{end+1,1} = 'Errors are assumed as additive Gaussian errors';
    if ar.config.fiterrors == -1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)==1)==0)
        General_string(end+2,1:2) = {'Chi2 value of the model is ', chi2Val.chi2_all};
        General_string(end+3,1:2)={'Compartments','size'};
    else
        General_string(end+2,1:2) = {'Log-likelihood value of the model is ', chi2Val.loglik_all};
        General_string(end+1,1:2) = {'Chi2 value of the model is ', chi2Val.chi2_res + chi2Val.chi2_prior + chi2Val.chi2_constr};
        General_string(end+2,1:2)={'Compartments','size'};
    end
    
    uses_ss = false;
    if(isfield(ar,'ss_conditions'))
        uses_ss = ar.ss_conditions;
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
    
    General_string{end+1,1} = 'Numerical steady-state equilibration';
    if(uses_ss==0)
        General_string{end,2} = 'no';
    else
        General_string{end,2} = 'yes';
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
        if (~exist('steadystate','var'))
            if ( isfield( ar.model(im), 'ss_condition' ) )
                steadystate = true;
            else
                steadystate = false;
            end
        end

        %loop over data to get condition values
        for jd= 1:length(ar.model(im).data)
            if(jd==1)
                pold_tmp = ar.model(im).data(jd).pold;
                continue;
            end
            not_inList = ~ismember(ar.model(im).data(jd).pold,pold_tmp);
            if(~isempty(not_inList) && sum(not_inList)>0)
                pold_tmp(end+1:end+sum(not_inList)) = ar.model(im).data(jd).pold(not_inList);
            end
        
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
        num_cols = 7;
        %Collect parameter trafos
        Conditions_string(end+2,1:num_cols) = {['Model file ' num2str(im)],'','exp condition','nTimePoints','Predictor','nDataPoints','chi2 value'}; 
        for jd= 1:length(ar.model(im).data)
           jc = ar.model(im).data(jd).cLink;   
           simulated_ss = 0;
            if ( steadystate )
                if ( ~isempty( ar.model(im).condition(jc).ssLink ) )
                    warning( 'Using simulated steady state values as initial condition (non-parametric)' );
                    x_ss = ar.model(im).ss_condition(ar.model(im).condition(jc).ssLink).xFineSimu(end,:);
                    simulated_states = ar.model(im).ss_condition(ar.model(im).condition(jc).ssLink).ssStates;
                    simulated_ss = 1;
                end
            end          
           
           Conditions_string{end+1,2} = ['Data file ' num2str(jd)];
           Conditions_string{end,3} = jc;
           Conditions_string{end,4} = length(ar.model(im).data(jd).tExp);
           if(strcmp(ar.model(im).data(jd).t,'t'))
               Conditions_string{end,5} = 'time';
           else
               Conditions_string{end,5} = ar.model(im).data(jd).t;
           end
           Conditions_string{end,6} = arnansum(ar.model(im).data(jd).ndata);
           Conditions_string{end,7} = arnansum(ar.model(im).data(jd).chi2);

           %append parameter trafos in each data struct
           [which_tmp,which_pold] = ismember(pold_tmp,ar.model(im).data(jd).pold);           
           par_trafo(which_tmp,jd) = regexprep(regexprep(ar.model(im).data(jd).fp(which_pold(which_pold~=0)),'^(',''),')$','');
           
           %Change values if steady-state equilibration was used
           if(simulated_ss)
               xName_tmp = regexprep(ar.model(im).x,'^.','init_$0');
               [which_tmp,which_ssState] = ismember(pold_tmp,xName_tmp(simulated_states==1));           
               par_trafo(which_tmp,jd) = cellfun(@num2str,num2cell(x_ss(which_ssState(which_ssState~=0))), 'UniformOutput', false);          
           end
           
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
   
    
    % Write out RAW ODEs and observation functions, to see parameter trafos and dependencies
    rawODE_string = {'ODE equations'};
    for im = 1:length(ar.model)
        rawODE_string{end+2,1} = ['Model ' num2str(im)];
        tmp_ode = ar.model(im).fx;
        for i = 1:length(ar.model(im).x)
            rawODE_string{end+1,1} = ['d' ar.model(im).x{i} '/dt'];                
            rawODE_string{end,2} = char(sym(tmp_ode{i}));                
        end
        
        if(isfield(ar.model(im),'yNames') && ~isempty(ar.model(im).yNames))
            rawODE_string{end+2,1} = 'Observables';
            rawODE_string(end+1,2:4) = {'scale','observation function','uncertainty parameter'};

            for i = 1:length(ar.model(im).yNames)                       
                if(sum(ar.model(im).logfitting(i))==0)
                    tmp_text = 'non-log';
                elseif(sum(ar.model(im).logfitting(i))==1)
                    tmp_text = 'log';
                end
               rawODE_string{end+1,1} = ar.model(im).yNames{i};
               rawODE_string{end,2} = tmp_text;

               tmp_ys = ar.model(im).fy{i};
               tmp_ystd = ar.model(im).fystd{i};

               if(strcmp(tmp_text,'log'))
                   rawODE_string{end,3} = ['log10(' tmp_ys ')'];
               else
                   rawODE_string{end,3} = tmp_ys;
               end
               rawODE_string{end,4} = tmp_ystd;

            end
        
            if(~isempty(ar.model(im).fz))
                found_z = 0;
                rawODE_string{end+2,1} = '';
                for jz = 1:length(ar.model(im).z)
                    regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(im).z{jz},ar.model(im).z{jz},ar.model(im).z{jz},ar.model(im).z{jz});
                    z_index = regexp(ar.model(im).data(id).fy,regPar);
                    check_empty = cellfun(@isempty,z_index);
                    if(sum(check_empty)<length(check_empty))
                        rawODE_string{end+1,1} = ar.model(im).z{jz};
                        rawODE_string{end,2} = ar.model(im).fz{jz};
                        found_z = found_z+1;
                    end
                end  
                if(found_z>0)
                    rawODE_string{end-found_z,1} = 'With definitions';
                end
            end
        end
    end   
    
    xlwrite(['./Benchmark_paper/General_info.xlsx'],rawODE_string,'Raw ODEs');   
    
    
%     return;
    %% Data specific csv files
    specialFunc = { ...
            {'smoothstep1',     '%s + (%s-%s) / (exp((%s-%s) / %s) + 1)', [4, 2, 4, 1, 3, 5], 'smoothstep1(t, level1, switch_time, level2, smoothness)' }, ...
            {'smoothstep2',     '%s + (%s-%s) / (exp((%s-%s) / %s) + 1) + (%s-%s) / (exp((%s-%s) / %s) + 1)', [6, 2, 4, 1, 3, 7, 4, 6, 1, 5, 7], 'smoothstep2( t, level1, switch_time1, level2, switch_time2, level3, smoothness )' }, ...        
            {'step1',           '%s + (%s-%s) * heaviside(%s-%s)', [2, 4, 2, 1, 3], 'step1(t, level1, switch_time, level2)'}, ...
            {'step2',           '%s + (%s-%s) * heaviside(%s-%s) + (%s-%s)*heaviside(%s-%s)', [2, 4, 2, 1, 3, 6, 4, 1, 5], 'step2(t, level1, switch_time1, level2, switch_time2, level3)' }, ...
            {'bolus',           '%s * (1 / sqrt( 2 * pi * %s^2 ) ) * exp(-(%s - %s)^2 / (2*%s^2))', [2, 4, 1, 3, 4], 'bolus(t, amount, time_point, duration)' }, ...
            {'hill_ka',         '%s^%s / (%s^%s + %s^%s)', [1, 3, 2, 3, 1, 3], 'hill_ka(conc, ka, n )' }, ...
            {'hill_kd',         '%s^%s / (%s + %s^%s)', [1, 3, 2, 1, 3], 'hill_kd(conc, kd, n )' }, ...
            {'isnonzero',       '(2*heaviside(%s))-1', 1, 'isnonzero(level)'}, ...
            {'max2',            '0.5*(%s+%s+abs(%s-%s))', [1,2,1,2], 'max2(a, b)'}, ...
            {'smooth1',         '(heaviside((%s - %s)/(%s - %s)) - (heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s))^2*((2*heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s) - 2*heaviside((%s - %s)/(%s - %s)) + 3)', [3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3], 'smooth1(t, start, end)' }, ...
            {'smooth2',         '-(heaviside((%s - %s)/(%s - %s)) - (heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s))^3*((heaviside((%s - %s)/(%s - %s)) - (heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s))*((6*heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s) - 6*heaviside((%s - %s)/(%s - %s)) + 15) - 10)', [3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3], 'smooth2(t, start, end)' }, ...
        };

    %ODE equations
    for im = 1:length(ar.model)
        if (~exist('steadystate','var'))
            if ( isfield( ar.model(im), 'ss_condition' ) )
                steadystate = true;
            else
                steadystate = false;
            end
        end
        ode_string = {};
        for id = 1:length(ar.model(im).data)
            file_name = ['./Benchmark_paper/Model/' 'model' num2str(im) '_data' num2str(id) '.xlsx'];

            if(~isempty(ar.model(im).x))
                jc = ar.model(im).data(id).cLink;
                %file_name = ['./Benchmark_paper/Model/' 'model' num2str(im) '_condition' num2str(jc) '.xlsx'];
                tmp_ode = ar.model(im).fx;
                fx = sym(tmp_ode);
                fus = regexprep(regexprep(ar.model(im).data(id).fu,'^.','($0'),'.$','$0)');
                fx = subs(fx, ar.model(im).u, fus');  
                
                fzs = regexprep(regexprep(ar.model(im).fz,'^.','($0'),'.$','$0)');
                fx = subs(fx, ar.model(im).z, fzs');         
                
                fx = subs(fx, ar.model(im).condition(jc).pold, ar.model(im).condition(jc).fp');      
                
                for i = 1:length(ar.model(im).x)
                    ode_string{i+2,1} = ['d' ar.model(im).x{i} '/dt'];
                    fx(i) = replaceFunctions( char(fx(i)), specialFunc, 0 );                
                    ode_string{i+2,2} = char(fx(i));                
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

            for i = 1:length(ar.model(im).data(id).y)            

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
               obs_string{end+1,1} = ar.model(im).data(id).y{i};
               obs_string{end,2} = tmp_text;
               
               tmp_ys = sym(ar.model(im).data(id).fy{i});
               tmp_ystd = sym(ar.model(im).data(id).fystd{i});
               
               fus = regexprep(regexprep(ar.model(im).data(id).fu,'^.','($0'),'.$','$0)');
               tmp_ys = subs(tmp_ys,ar.model(im).u,fus');
               tmp_ys = subs(tmp_ys,ar.model(im).data(id).pold,ar.model(im).data(id).fp');
               
               tmp_ystd = subs(tmp_ystd,ar.model(im).data(id).pold,ar.model(im).data(id).fp');
               if(strcmp(tmp_text,'log'))
                   obs_string{end,3} = ['log10(' char(tmp_ys) ')'];
               else
                   obs_string{end,3} = char(tmp_ys);
               end
               obs_string{end,4} = A_tmp;
               obs_string{end,5} = char(tmp_ystd);

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
                        tmp_zf = sym(ar.model(im).fz{jz});
                        tmp_zf = subs(tmp_zf,ar.model(im).data(id).pold,ar.model(im).data(id).fp');
               
                        obs_string{end,2} = char(tmp_zf);
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
            init_string{end+1,1} = 'Integration start';
            init_string{end,2} = ar.model(im).condition(jc).tstart;
            tmp = ar.model(im).px0;
            simulated_ss = 0;
            if ( steadystate )
                if ( ~isempty( ar.model(im).condition(jc).ssLink ) )
                    warning( 'Using simulated steady state values as initial condition (non-parametric)' );
                    x_ss = ar.model(im).ss_condition(ar.model(im).condition(jc).ssLink).xFineSimu(end,:);
                    simulated_ss = 1;
                end
            end
            
                
            for i = 1:length(tmp)
               init_string{end+1,1} = tmp{i};
               init_tmp = sym(ar.model(im).condition(jc).fp{ismember(ar.model(im).condition(jc).pold,tmp{i})});
               fzs = ar.model(im).fz; %regexprep(regexprep(ar.model(im).fz,'^.','($0'),'.$','$0)');
               init_string{end,2} = subs(init_tmp, ar.model(im).z, fzs');
               if(~simulated_ss)
                  init_string{end,2} = char(init_tmp);  
               else
                   jx = strcmp(ar.model(im).x,strrep(tmp{i},'init_',''));
                   init_string{end,2} = x_ss(jx); 
               end
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

