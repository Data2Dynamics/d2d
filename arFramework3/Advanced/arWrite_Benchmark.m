function arWrite_Benchmark()
global ar
    if(~exist('./Benchmark_paper', 'dir'))
        mkdir('./Benchmark_paper')
    end
    
    %% General Infos
    
    A = {ar.model(1).name; ''; 'FACTS'};
    tmp = sprintf('The model contains %i data points %i free parameters and %i experimental conditions',ar.ndata,sum(ar.qFit==1),length(ar.model(1).condition));
    A = [A; tmp; ''; 'Observables'];
    if(isfield(ar.model(1),'logfitting'))
        if(sum(ar.model(1).logfitting)==0)
            tmp_text = 'non-log';
        elseif(sum(ar.model(1).logfitting)==length(ar.model(1).logfitting))
            tmp_text = 'log';
        else
            tmp_text = 'CHECK';
        end
        tmp = sprintf('%i observables are specified on %s scale',length(ar.model(1).fy),tmp_text);
        A = [A; tmp];
    end
    fid = fopen('./Benchmark_paper/general_info.csv','w');
    for i = 1:length(A)
        fprintf(fid, '%s\n',A{i});
    end
    fclose(fid)
    
    %% ODE equations
    
    A = {'ODE equations'; ''};
    
    for i = 1:length(ar.model(1).x)
        A{i+2,1} = ['d' ar.model(1).x{i} '/dt'];
        tmp_ode = '';
        for j = 1:size(ar.model(1).N,2)
            if(ar.model(1).N(i,j)==1)
                tmp_ode = [tmp_ode ' + ' ar.model(1).fv{j}];
            elseif(ar.model(1).N(i,j) == -1)
                tmp_ode = [tmp_ode ' - ' ar.model(1).fv{j}];
            elseif(ar.model(1).N(i,j) ~= 0)
                tmp_ode = [tmp_ode ' ' ar.model(1).N(i,j) ' * ' ar.model(1).fv{j}];
            end
        end
        
        for j = 1:length(ar.model(1).x)
            regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(1).x{j},ar.model(1).x{j},ar.model(1).x{j},ar.model(1).x{j});
            tmp_ode = regexprep(tmp_ode,regPar, ['[' ar.model(1).x{j} ']']);
        end
        A{i+2,2} = tmp_ode;
    end 
    
    fid = fopen('./Benchmark_paper/odes.csv','w');
    for i = 1:size(A,1)
        fprintf(fid, '%s, %s \n',A{i,1},A{i,2});
    end
    fclose(fid)
    
    %% Parameter values
    
    A = {'Parameter values';'parameter'};
    A(2,2:5) = {'value','lower boundary','upper boundary','analysis at log-scale'};
    
    for i = 1:length(ar.p)
        A{i+2,1} = ar.pLabel{i}; 
        A{i+2,2} = num2str(ar.p(i));
        A{i+2,3} = num2str(ar.lb(i));
        A{i+2,4} = num2str(ar.ub(i));
        A{i+2,5} = num2str(ar.qLog10(i));
    end
    
    fid = fopen('./Benchmark_paper/parameters.csv','w');
    for i = 1:size(A,1)
        fprintf(fid, '%s, %s, %s, %s, %s \n',A{i,1},A{i,2},A{i,3},A{i,4},A{i,5});
    end
    fclose(fid)
    
    %% Conditions
    A = {'Conditions';''};   
    tmp_par = {};
    for jc= 1:length(ar.model(1).condition)
       A{jc+2,1} = ['Condition' num2str(jc)]; 
        tmp_id = ar.model(1).condition(jc).dLink(1);
        
        for jd=1:length(ar.model(1).data(tmp_id).condition)
            if(sum(strcmp(ar.model(1).data(tmp_id).condition(jd).parameter,tmp_par))==0)
                tmp_par = [tmp_par ar.model(1).data(tmp_id).condition(jd).parameter];
            end
        end
    end
    A(2,2:length(tmp_par)+1) = tmp_par;
    for jc= 1:length(ar.model(1).condition)
        tmp_id = ar.model(1).condition(jc).dLink(1);
        for jd=1:length(ar.model(1).data(tmp_id).condition)
            which_cond = find(ismember(tmp_par,ar.model(1).data(tmp_id).condition(jd).parameter));
            A{jc+2,which_cond+1} = num2str(ar.model(1).data(tmp_id).condition(jd).value);
        end

    end       
   
    fid = fopen('./Benchmark_paper/conditions.csv','w');
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            fprintf(fid, '%s,',A{i,j});
            
            if(j==size(A,2))
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid)
    
    %% Observations
    
    A = {'Observables';''};
    
    A(2,2:5) = {'scale','measurement','observation function','uncertainties'};
    for id = 1:length(ar.model(1).data)
        for i = 1:length(ar.model(1).data(id).yNames)
            found = 0;
            for itest = 2:size(A,1)
                if(strcmp(ar.model(1).data(id).yNames{i},A{itest,1}))
                    found = 1;
                end
            end
            if(found)
                continue;
            end
            if(sum(ar.model(1).data(id).logfitting(i))==0)
                tmp_text = 'non-log';
            elseif(sum(ar.model(1).data(id).logfitting(i))==1)
                tmp_text = 'log';
            end
           A{end+1,1} = ar.model(1).data(id).yNames{i};
           A{end,2} = tmp_text;
           A{end,4} = ar.model(1).data(id).fy{i};
           if(sum(isnan(ar.model(1).data(id).yExpStd(:,i)))==length(ar.model(1).data(id).tExp))
               A{end,5} = 'fitted';
           else
               A{end,5} = 'measured';
           end
        end
    end
    
    A{end+2,1} = 'With definitions';
    for i = 1:length(ar.model(1).fz)
       A{end+1,1} = ar.model(1).z{i};
       A{end,2} = ar.model(1).fz{i};
    end
    
    A{end+2,1} = 'Error model:';
    A{end+1,1} = 'Fill out manually!!';
    
    fid = fopen('./Benchmark_paper/observables.csv','w');
    for i = 1:size(A,1)
        fprintf(fid, '%s, %s, %s, %s, %s \n',A{i,1},A{i,2},A{i,3},A{i,4},A{i,5});
    end
    fclose(fid)
    
    %% Initials
    
    A = {'Initial values';''};
    tmp = ~cellfun(@isempty,strfind(ar.model(1).condition(1).pold,'init_'));
    
    for i = find(tmp)
       A{end+1,1} = ar.model(1).condition(1).pold{i};
       A{end,2} = strrep(strrep(ar.model(1).condition(1).fp{i},'(',''),')','');
    end
    
    fid = fopen('./Benchmark_paper/initials.csv','w');
    for i = 1:size(A,1)
        fprintf(fid, '%s, %s \n',A{i,1},A{i,2});
    end
    fclose(fid)
    
    
    %% Constants
    
    A = {'Constant values';''};
    [tmp,ia] = setdiff(ar.model(1).condition(1).pold,ar.model(1).condition(1).fp);
    ia(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
    tmp(~cellfun(@isempty,strfind(tmp,'init_'))) = [];
    for i = 1:length(tmp)
       A{end+1,1} = tmp{i};
       A{end,2} = strrep(strrep(ar.model(1).condition(1).fp{ia(i)},'(',''),')','');
    end
    
    fid = fopen('./Benchmark_paper/constants.csv','w');
    for i = 1:size(A,1)
        fprintf(fid, '%s, %s \n',A{i,1},A{i,2});
    end
    fclose(fid)
end