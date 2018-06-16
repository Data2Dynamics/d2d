function arWrite_Benchmark(im)
global ar
    if(~exist('./Benchmark_paper', 'dir'))
        mkdir('./Benchmark_paper')
    end
    if(~exist('./Benchmark_paper/Model', 'dir'))
        mkdir('./Benchmark_paper/Model')
    end
    if(~exist('im','var') || isempty(im))
        im = 1;
    end
    %% General csv sheets
    
    %General Infos
    
    General_string = {ar.model(im).name; ''; 'FACTS'};
    tmp = sprintf('The model contains %i data points %i free parameters and %i experimental conditions',ar.ndata,sum(ar.qFit==1),length(ar.model(im).condition));
    General_string = [General_string; tmp];
    General_string = [General_string; 'Errors are assumend as additive Gaussian errors'];
    General_string(end+2,1:2) = {'Chi2 value of the model is ', ar.chi2};
    General_string(end+3,1:2)={'Compartments','size'};
    for i = 1:length(ar.model(im).c)
        General_string = [General_string; {ar.model(im).c{i}, str2num(ar.model(im).pc{i})}];
    end
    
    General_string{end+2,1} = 'Special features of the model';
    General_string{end+1,1} = 'Splines';
    if(contains(ar.model(im).fu,'spline'))
        General_string{end,2} = 'yes';
    else
        General_string{end,2} = 'no';
    end
    
    General_string{end+1,1} = 'Estimated error parameters';
    est_errors = 'no';
    log_fit = 'no';
    for id = 1:length(ar.model(im).data)
        if(sum(sum(isnan(ar.model(im).data(id).yExpStd)))>sum(sum(isnan(ar.model(im).data(id).yExp))) && sum(sum(isnan(ar.model(im).data(id).yExp)))<(length(ar.model(im).data(id).tExp)*size(ar.model(im).data(id).yExp,2)))
           est_errors = 'yes';     
        end
        if(sum(ar.model(im).data(id).logfitting)>0)
           log_fit = 'yes'; 
        end
    end
    General_string{end,2} = est_errors;
    
    if(strcmp(est_errors,'yes'))
       General_string{end,3} = 'Check sheet Parameters_noErrorModel for values with fixed errors'; 
    end
    
    General_string{end+1,1} = 'log-scale of observations';
    General_string{end,2} = log_fit;
    
    General_string{end+1,1} = 'Events';
    General_string{end,2} = 'Fill out manually?';
    
    General_string{end+1,1} = 'mixed-effects';
    General_string{end,2} = 'No';   
        
    
    xlwrite('./Benchmark_paper/General_info.xlsx',General_string,'General Info');
    
     % Parameter values
    
    General_string = {'Parameter values';'parameter'};
    General_string(2,2:5) = {'value','lower boundary','upper boundary','analysis at log-scale'};
    
    for i = 1:length(ar.p)
        General_string{i+2,1} = ar.pLabel{i}; 
        General_string{i+2,2} = ar.p(i);
        General_string{i+2,3} = ar.lb(i);
        General_string{i+2,4} = ar.ub(i);
        General_string{i+2,5} = ar.qLog10(i);
    end
    
    xlwrite('./Benchmark_paper/General_info.xlsx',General_string,'Parameters');
    
    % Conditions
    Conditions_string = {'Conditions','';'','exp condition'};   
%     tmp_par = {};
    max_fp = 1;
    max_id = 1;
    for jc = 1:length(ar.model(im).condition)
        if(length(ar.model(im).condition(jc).fp)>max_fp)
           max_fp = length(ar.model(im).condition(jc).fp);
           max_id = jc;
        end
    end
    par_trafo = cell(max_fp,length(ar.model(im).data));
    %Collect parameter trafos
    for jd= 1:length(ar.model(im).data)
        jc = ar.model(im).data(jd).cLink;
       Conditions_string{jd+2,1} = ['Data file ' num2str(jd)];
       Conditions_string{jd+2,2} = num2str(jc);
       %append parameter trafos in each data struct
       par_trafo(1:length(ar.model(im).condition(jc).fp),jd) = ar.model(im).condition(jc).fp;
    end
    
    %get differences in parameter trafos
    for i = 1:size(par_trafo,1)
        ids_empty = cellfun(@isempty,par_trafo(i,:));
        if(any(ids_empty))
            for ie = find(ids_empty==1)
                par_trafo{i,ie} = '';
            end
        end
        num(i) = length(unique(par_trafo(i,:)));
        istrafo(i) = num(i)==1 & ~strcmp(ar.model(1).condition(max_id).fp{i},ar.model(1).condition(max_id).pold{i});
    end   
    tmp_par = ar.model(1).condition(max_id).pold(num>1);
    
    %Go through differing parameter trafos (num>1)
    Conditions_string(2,3:length(tmp_par)+2) = tmp_par;
    Conditions_string(3:length(ar.model(im).data)+2,3:length(tmp_par)+2) = par_trafo(num>1,:)';%strrep(strrep(par_trafo(num>1,:)','(',''),')','');
    
    Conditions_string(2,end+1:end+3) = {'nTimePoints','nDataPoints','chi2 value'};
    for jd= 1:length(ar.model(im).data)
        Conditions_string{jd+2,end-2} = length(ar.model(im).data(jd).tExp);
        Conditions_string{jd+2,end-1} = arnansum(ar.model(im).data(jd).ndata);
        Conditions_string{jd+2,end} = arnansum(ar.model(im).data(jd).chi2);
    end
    %This part writes out the constants (num==1 means same parameter trafo
    %in every condition)
    Conditions_string(end+2,1) = {'General transformations'};
    
    Conditions_string(end+1:size(Conditions_string,1)+sum(istrafo),2) = ar.model(im).condition(jc).pold(istrafo)';   
    Conditions_string(size(Conditions_string,1)-sum(istrafo)+1:size(Conditions_string,1),3) = ar.model(im).condition(jc).fp(istrafo);%strrep(strrep(ar.model(im).condition(jc).fp(istrafo),'(',''),')','');        
   
    xlwrite('./Benchmark_paper/General_info.xlsx',Conditions_string,'Experimental conditions');
   
    
    %% Data specific csv files
    
    %ODE equations
    
    for id = 1:length(ar.model(im).data)
        file_name = ['./Benchmark_paper/Model/' ar.model(im).name '_model_data' num2str(id) '.xlsx'];
    
        if(~isempty(ar.model(im).x))
            tmp_ode = ar.model(im).fx;
            for j = 1:length(ar.model(im).u)
                regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(im).u{j},ar.model(im).u{j},ar.model(im).u{j},ar.model(im).u{j});
                tmp_ode = regexprep(tmp_ode,regPar, ar.model(im).fu{j});    
            end
            
            jc = ar.model(im).data(id).cLink;
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

        obs_string(2,2:6) = {'scale','observation function','uncertainties','error model','data max normalized to 1?'};
    
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
           
           obs_string{end,6} = ar.model(im).data(id).normalize(i);
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
