function arWrite_Benchmark(im)
global ar
    if(~exist('./Benchmark_paper', 'dir'))
        mkdir('./Benchmark_paper')
    end
    if(~exist('im','var') || isempty(im))
        im = 1;
    end
    %% General Infos
    
    A = {ar.model(im).name; ''; 'FACTS'};
    tmp = sprintf('The model contains %i data points %i free parameters and %i experimental conditions',ar.ndata,sum(ar.qFit==1),length(ar.model(im).condition));
    A = [A; tmp; ''; 'Observables'];
    if(isfield(ar.model(im),'logfitting'))
        if(sum(ar.model(im).logfitting)==0)
            tmp_text = 'non-log';
        elseif(sum(ar.model(im).logfitting)==length(ar.model(im).logfitting))
            tmp_text = 'log';
        else
            tmp_text = 'CHECK';
        end
        tmp = sprintf('%i observables are specified on %s scale',length(ar.model(im).fy),tmp_text);
        A = [A; tmp];
    end
    fid = fopen(['./Benchmark_paper/general_info_m' num2str(im) '.csv'],'w');
    for i = 1:length(A)
        fprintf(fid, '%s\n',A{i});
    end
    fclose(fid)
    
    %% ODE equations
    if(~isempty(ar.model(im).x))
        A = {'ODE equations'; ''};

        for i = 1:length(ar.model(im).x)
            A{i+2,1} = ['d' ar.model(im).x{i} '/dt'];
            A{i+2,2} = ar.model(im).fx_par{i};
    %         tmp_ode = '';
    %         for j = 1:size(ar.model(im).N,2)
    %             if(ar.model(im).N(i,j)==1)
    %                 tmp_ode = [tmp_ode ' + ' ar.model(im).fv{j}];
    %             elseif(ar.model(im).N(i,j) == -1)
    %                 tmp_ode = [tmp_ode ' - ' ar.model(im).fv{j}];
    %             elseif(ar.model(im).N(i,j) ~= 0)
    %                 tmp_ode = [tmp_ode ' ' ar.model(im).N(i,j) ' * ' ar.model(im).fv{j}];
    %             end
    %         end
    %         
    %         for j = 1:length(ar.model(im).x)
    %             regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.model(im).x{j},ar.model(im).x{j},ar.model(im).x{j},ar.model(im).x{j});
    %             tmp_ode = regexprep(tmp_ode,regPar, ['[' ar.model(im).x{j} ']']);
    %         end
    %         A{i+2,2} = tmp_ode;
        end 
              
        % Inputs    
        if(~isempty(ar.model(im).fu))

            A(length(ar.model(im).x)+3:length(ar.model(im).x)+5,1:2) = {'', '';'The following', ' parameters are obtained from input functions';'Parameter', 'Input function'};

            for i=1:length(ar.model(im).u)
                A{end+1,1} = ar.model(im).u{i};
                A{end,2} = ar.model(im).fu{i};
            end
        end

        fid = fopen(['./Benchmark_paper/odes_m' num2str(im) '.csv'],'w');
        for i = 1:size(A,1)
            fprintf(fid, '%s; %s \n',A{i,1},A{i,2});
        end
        fclose(fid)
    end
    
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
    
    fid = fopen(['./Benchmark_paper/parameters_m' num2str(im) '.csv'],'w');
    for i = 1:size(A,1)
        fprintf(fid, '%s; %s; %s; %s; %s \n',A{i,1},A{i,2},A{i,3},A{i,4},A{i,5});
    end
    fclose(fid)
    
    %% Conditions
    A = {'Conditions';''};   
%     tmp_par = {};
    hae = {};
    for jc= 1:length(ar.model(im).condition)
       A{jc+2,1} = ['Condition' num2str(jc)];
       hae = [hae ar.model(im).condition(jc).fp];
%         tmp_id = ar.model(im).condition(jc).dLink(1);        
%         for jd=1:length(ar.model(im).data(tmp_id).condition)
%             if(sum(strcmp(ar.model(im).data(tmp_id).condition(jd).parameter,tmp_par))==0)
%                 tmp_par = [tmp_par ar.model(im).data(tmp_id).condition(jd).parameter];
%             end
%         end
    end
    
    for i = 1:size(hae,1)
       num(i) = length(unique(hae(i,:))); 
    end    
    tmp_par = ar.model(1).condition(1).pold(num>1);
    
    A(2,2:length(tmp_par)+1) = tmp_par;
    for icond = find(num>1)
        A(3:length(ar.model(im).condition)+2,2:length(tmp_par)+1) = strrep(strrep(hae(num>1,:)','(',''),')','');
    end
%     for jc= 1:length(ar.model(im).condition)
%         tmp_id = ar.model(im).condition(jc).dLink(1);
%         for jd=1:length(ar.model(im).data(tmp_id).condition)
%             which_cond = find(ismember(tmp_par,ar.model(im).data(tmp_id).condition(jd).parameter));
%             A{jc+2,which_cond+1} = num2str(ar.model(im).data(tmp_id).condition(jd).value);
%         end
%     end       
   
    fid = fopen(['./Benchmark_paper/conditions_m' num2str(im) '.csv'],'w');
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            fprintf(fid, '%s;',A{i,j});
            
            if(j==size(A,2))
                fprintf(fid,'\n');
            end
        end
    end
    fclose(fid)
    
    % Constants
    A_first = {'Constant values';'Check if constants are representative for model!'};
    tmp_nr = find(~ismember(ar.model(1).condition(1).pold, hae(:,1)') & num==1 & cellfun(@isempty,strfind(ar.model(1).condition(1).pold,'init_')));
    A = ar.model(1).condition(1).pold(tmp_nr)';
    A(:,2) = strrep(strrep(hae(tmp_nr,1),'(',''),')','');
    
    fid = fopen(['./Benchmark_paper/constants_m' num2str(im) '.csv'],'w');
    fprintf(fid,'%s; \n %s; \n',A_first{1},A_first{2});
    for i = 1:size(A,1)
        fprintf(fid, '%s; %s \n',A{i,1},A{i,2});
    end
    fclose(fid)
    
    %% Observations
    
    A = {'Observables';''};
    
    A(2,2:5) = {'scale','measurement','observation function','uncertainties'};
    for id = 1:length(ar.model(im).data)
        for i = 1:length(ar.model(im).data(id).yNames)
            found = 0;
            if(sum(strcmp(A(2:end,1),ar.model(im).data(id).yNames{i}))>0)
                found = 1;
                tmp = A{strcmp(A(2:end,1),ar.model(im).data(id).yNames{i}),5};
            end
            
            for itest = 2:size(A,1)
                if(strcmp(ar.model(im).data(id).yNames{i},A{itest,1}))
                    found = 1;
                end
            end
            if(sum(isnan(ar.model(im).data(id).yExpStd(:,i)))==length(ar.model(im).data(id).tExp) && sum(isnan(ar.model(im).data(id).yExp(:,i)))<length(ar.model(im).data(id).tExp))
               A_tmp = 'fitted';
            elseif(sum(isnan(ar.model(im).data(id).yExpStd(:,i)))<length(ar.model(im).data(id).tExp) && sum(isnan(ar.model(im).data(id).yExp(:,i)))<length(ar.model(im).data(id).tExp))
               A_tmp = 'measured';
            else
               continue;
            end
            if(found && (strcmp(A_tmp,'fitted') && strcmp(tmp,'measured')) || (strcmp(A_tmp,'measured') && strcmp(tmp,'fitted')))
               found = 0;
               A_tmp = 'mixed';
            end
            if(found)
                continue;
            end
            if(sum(ar.model(im).data(id).logfitting(i))==0)
                tmp_text = 'non-log';
            elseif(sum(ar.model(im).data(id).logfitting(i))==1)
                tmp_text = 'log';
            end
           A{end+1,1} = ar.model(im).data(id).yNames{i};
           A{end,2} = tmp_text;
           A{end,4} = ar.model(im).data(id).fy{i};
           A{end,5} = A_tmp;
        end
    end
    
    if(~isempty(ar.model(im).fz))
        A{end+2,1} = 'With definitions';
        for i = 1:length(ar.model(im).fz)
           A{end+1,1} = ar.model(im).z{i};
           A{end,2} = ar.model(im).fz{i};
        end
    end
    
    A{end+2,1} = 'Error model:';
    A{end+1,1} = 'Fill out manually!!';
    
    fid = fopen(['./Benchmark_paper/observables_m' num2str(im) '.csv'],'w');
    for i = 1:size(A,1)
        fprintf(fid, '%s; %s; %s; %s; %s \n',A{i,1},A{i,2},A{i,3},A{i,4},A{i,5});
    end
    fclose(fid)
    
    %% Initials
    
    A = {'Initial values';''};
    tmp = ~cellfun(@isempty,strfind(ar.model(im).condition(1).pold,'init_'));
    
    for i = find(tmp)
       A{end+1,1} = ar.model(im).condition(1).pold{i};
       A{end,2} = strrep(strrep(ar.model(im).condition(1).fp{i},'(',''),')','');
    end
    
    fid = fopen(['./Benchmark_paper/initials_m' num2str(im) '.csv'],'w');
    for i = 1:size(A,1)
        fprintf(fid, '%s; %s \n',A{i,1},A{i,2});
    end
    fclose(fid)
    
    
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