%% AUTOMATIC REPARAMETRIZATION
% 1.- Load data.
% 2.- Compute the necessary number of reparametrizations.
% 3.- First study of observability and symmetries.
% 4.- Reparametrization loop:
%   4.1.- Type of generator and symmetry.
%   4.2.- Choose the parameter to remove.
%   4.3.- Reparametrization of the model:
%       4.3.1.- Isolate exp(epsilon)|| epsilon.
%       4.3.2.- Replace it in the vector of new variables.
%       4.3.3.- Obtain the reparametrized model
%   4.4.- Redefinition of states and parameters. Save the model and
%         compute the observability and symmetry study.
%
%

function AutoRepar

%clearvars ar
%clc

global ar


% LOAD DATA
%   Path for data  
% parentpath = cd(cd('..'));
% if 0==isfolder(strcat(parentpath,filesep,'functions'))
%     parentpath=strcat(pwd,filesep,'functions');
%     cd(parentpath) 
%     fprintf('Your folder has been changed. \n')
% end
parentpath = ar.ia.paths.functions;
modelspath = ar.ia.paths.models;
resultspath= ar.ia.paths.results;
cd(parentpath)

% addpath(parentpath);
% addpath(modelspath);
% addpath(resultspath);



%% STRIKE-GOLDD
[modelname,paths,opts,submodels,prev_ident_pars] = options();
if (opts.use_existing_results==0)
    STRIKE_GOLDD(modelname,1);
    resultsname = sprintf('id_results_%s_%s',modelname,date);
    load(resultsname);
else
    load(opts.results_file);
    [modelname,paths,opts,submodels,prev_ident_pars] = options();
end

%% 
syms epsilon epsilon_t 
assume(epsilon,'real');
if exist('w')~=1
   w=[]; 
elseif exist('u')~=1
   u=[];
end
%% NUMBER OF REPARAMETRIZATIONS
num_repar = length(xaug)-rango;
if num_repar==0
    fprintf('\n-----------------------------------------\n');
    fprintf('The model is already FISPO.\n There is no need for reparameterization.\n')
    fprintf('-----------------------------------------\n\n');

    return
else
    fprintf('\n--------------------------------------------------\n');
    fprintf('>>> The %s model requires %d transformation(s).\n', modelname,num_repar);
end

%% SYMMETRY STUDY
fprintf('Searching for symmetries of %s...\n', modelname);
[transf,nuevas_variables,allVar,z_v]=Lie_Symmetry();

% Verify that the transformation vector is not empty
if (isempty(transf)==1)
   return 
end

%%  REPARAMETRIZATION LOOP
for it=1:num_repar
    fprintf('\n--------------------------------------------------');
    fprintf('\n>>> Reparameterization #%d:\n', it);
    l_allVar    =   length(allVar);
    transf_new  =   transf;
    nv_new      =   simplify(nuevas_variables);
    [c,r]       =   size(transf_new);
    A           =   zeros(c,r);
    cnbs        =   []; % New variables vector 
    cols        =   []; % Row position of transformations 
    dim         =   []; % Number of transformations for each infinitesimal generator
    txw         =   length(x)+length(w);
    tics        =   length(p)+txw;

    %% SYMMETRY TYPE AND PARAMETERS 
    for j=1:c
       [~,col]=find(transf_new(j,1:end));
       [~,c1]=find(col>txw & col<=tics);
       col=col(c1);
       cnb=nv_new(j,col);
       l=length(cnb);
       [num,~]=numden(cnb);
       num=expand(num);
       ind=[];
       for i=1:l
           temp1=children(collect(num(i),epsilon));
           l_temp1=length(temp1);
           if l_temp1 <3
               for k=1:l_temp1
                   [cc,tt]=coeffs(temp1(k),allVar);
                    if (cc==exp(epsilon) || cc==exp(-epsilon))
                        A(j,col(i))=1;  %   Scaling
                    else
                        A(j,col(i))=2;  %   Other type of elemenatry transf.
                    end           
               end
           else
              ind=[ind,i];
           end
       end
       col(ind)=[];
       cnb(ind)=[];
       %% DELETE INFINITESIMAL GENERATOR OF IDENTIFIABLE PARAMATERS
       if (exist('p_id')==1 && isempty(p_id)==0 && isempty(col)==0)
           [ind_p_id,~]=find(allVar==p_id);
           [~,ind_p_id_col]=find(col==ind_p_id);
           if (isempty(ind_p_id_col)==0)
                col(ind_p_id_col)=[];
           end
       end
       if (isempty(col)==0)
           cnbs=[cnbs,cnb];
           cols=[cols,col];
       end
       fprintf('   Generator #%d has %d possible reparameterization(s) \n',j,length(col));
       dim=[dim;length(col)];
    end
    [row1,~]=find(A==1);    %   Positions with scaling transf.
    [row2,~]=find(A==2);    %   Positions with other type of elementary transf.
    [val,~]=intersect(row1,row2);   %   Both
    l_val=length(val);

    %% CHOOSE POSSIBLE PARAMETERS TO REMOVE
    sum_num_par=[];
    for j=1:c
       fprintf('\n   For the infinitesimal generator #%d the following parameters can be removed:',j);
       for i=1:dim(j)
                if j~=1
                    fprintf('\n   (%d): %s  ',i,allVar(cols(i+sum(dim(1:j-1)))))
                    fprintf('\n     --> Affected variables: ')
                    [~,jj,~]=find(transf(j,:));
%                     disp(transpose(allVar(jj)))
                    fprintf('%s ',char(allVar(jj)))
                else
                    fprintf('\n   (%d): %s ',i,allVar(cols(i)))
                    fprintf('\n      --> Affected variables: ')
                    [~,jj,~]=find(transf(j,:));
%                     disp(transpose(allVar(jj)))
                    fprintf('%s ',char(allVar(jj)))
                end 
       end
        sum_num_par=[sum_num_par;dim(j)];
        fprintf('\n')
    end
    %   If no generators are available, then return
    if (sum_num_par==0)
        fprintf('No reparameterization is possible. Only states or known inputs can be removed. \n')
        return
    end
    %   Choose infinitesimal generator
    if length(sum_num_par)~=1
        prompt   =   '\n   Enter the generator that you want to use:\n';
        num_gen  =   input(prompt);
        %   Verification
        while ((num_gen==0)|| (num_gen>c) || (sum_num_par(num_gen)==0))
            fprintf('The choosen generator is not correct. Please choose another one. \n')
            prompt   =   '   Enter the generator that you want to use:\n';
            num_gen  =   input(prompt);
        end
    else
        num_gen=1;
    end
    if sum_num_par(num_gen)~=1
        %   Choose parameter
        prompt   =   '   Enter the number of the parameter to be removed (in order of output per screen):\n';
        num_par  =   input(prompt);
        %   Verification
        while ((num_par==0) || sum_num_par(num_gen)<num_par)
            fprintf('   The chosen parameter is not correct. Please choose another one. \n')
            prompt   =   '   Enter the number of the parameter to be removed:\n';
            num_par  =   input(prompt);
        end
    else
        num_par=sum_num_par(num_gen);
    end
    %   Imprimir por pantalla el par�metro del generador que se elimina
    fprintf('   You have chosen the parameter ')
    if num_gen ~=1
        fprintf('%s',allVar(cols(num_par+sum(dim(1:num_gen-1)))))
    else
        fprintf('%s',allVar(cols(num_par)))
    end 
    fprintf(' from the infinitesimal generator:\n')
    disp(transf(num_gen,:))

    %% REPARAMETRIZATION OF THE MODEL
    % 1.- Classify the transformation type
    % 2.- Isolate exp(epsilon) or epsilon
    % 3.- Reformulate the model

    if num_gen~=1
        cl  =   find(allVar(cols(num_par+sum(dim(1:num_gen-1))))==allVar);
    else
        cl  =   find(allVar(cols(num_par))==allVar);
    end
    if ( A(num_gen,cl)==1 )
        % Scaling
        if num_gen~=1
            eq1=nv_new(num_gen,cl)==1;
            epsilon_t=isolate(eq1,exp(epsilon));
            aux3=allVar;
            aux3(cols(num_par+sum(dim(1:num_gen-1))))=1;
        else
            eq1=nv_new(num_gen,cl)==1;
            epsilon_t=isolate(eq1,exp(epsilon));
            aux3=allVar;
            aux3(cols(num_par))=1;
        end
        epsilon_a=epsilon_t;
        chh=children(epsilon_t);
        if iscell(chh)==1
            chh=[chh{:}];
        end
        epsilon_t=chh(1,1);
        leq1=simplify(log(epsilon_t),'IgnoreAnalyticConstraints',true);
        ceq1=coeffs(leq1,epsilon);
        aux1=subs(nv_new(num_gen,:),exp(-ceq1*epsilon), 1/chh(1,2));
        aux2=subs(aux1, exp(ceq1*epsilon),chh(1,2));
        %------------------------------------------------------------------
        [~,col_z_v_3]=find(z_v(num_gen,:)==3);
        [~,col_z_v_2]=find(z_v(num_gen,:)==2);
        if (opts.ode_n==0 && isempty(col_z_v_3)==0)
            [exp_e,exp_p] = expand_lie(opts.tmax,allVar(cl));
            aux2(col_z_v_3)=simplify(subs(aux2(col_z_v_3),exp_e,exp_p));
        elseif isempty(col_z_v_2)==0
            chh_2=children(isolate(epsilon_a,epsilon));
            if iscell(chh_2)==1 % For R2020b
                chh_2=[chh_2{:}];
            end
            aux2=subs(aux2,epsilon,chh_2(1,2));
        end
        %------------------------------------------------------------------
        aux2=simp(allVar,aux2,chh(1,2));
        %------------------------------------------------------------------
        f_s=simplify(subs(f,allVar,aux3));
        h_s=simplify(subs(h,allVar,aux3));
        fprintf('>>> Transformed variables: \n')
        for ii=1:length(aux2)
            fprintf('%s <-- %s \n',allVar(ii),aux2(ii))
        end
        fprintf('>>> Equations of the reformulated model: \n')
        fprintf('f:\n')
        disp(f_s)
        fprintf('h:\n')
        disp(h_s)
    elseif ( A(num_gen,cl)==2 )
        % Other elementary transformations (excluding scaling)
        if num_gen~=1
            eq1=cnbs(num_par+sum(dim(1:num_gen-1)))==1;
            epsilon_t=isolate(eq1,epsilon);
            aux3=allVar;
            aux3(cols(num_par+sum(dim(1:num_gen-1))))=1;
        else
            eq1=cnbs(num_par)==1;
            epsilon_t=isolate(eq1,epsilon);
            aux3=allVar;
            aux3(cols(num_par))=1;
        end
        chh=children(epsilon_t);
        if iscell(chh)==1
            chh=[chh{:}];
        end
        epsilon_t=chh(1,2);
        aux1=subs(nv_new(num_gen,:),epsilon,epsilon_t);
        %------------------------------------------------------------------
        [~,col_z_v_3]=find(z_v(num_gen,:)==3);
        [~,col_z_v_2]=find(z_v(num_gen,:)==2);
        if (opts.ode_n==0 && isempty(col_z_v_3)==0)
            [exp_e,exp_p] = expand_lie(opts.tmax,allVar(cl));
            aux2(col_z_v_3)=simplify(subs(aux2(col_z_v_3),exp_e,exp_p));
        elseif isempty(col_z_v_2)==0
            chh_2=children(isolate(epsilon_a,epsilon));
            if iscell(chh_2)==1 % For R2020b
                chh_2=[chh_2{:}];
            end
            aux2=subs(aux2,epsilon,chh_2(1,2));
        end
        %------------------------------------------------------------------
        aux1=simp(allVar,aux1,chh(1,2));
        %------------------------------------------------------------------
        f_s=simplify(subs(f,allVar,aux3));
        h_s=simplify(subs(h,allVar,aux3));
        fprintf('>>> Transformed variables: \n')
        for ii=1:length(aux1)
            fprintf('%s <-- %s \n',allVar(ii),aux1(ii))
        end
        fprintf('>>> Reformulated model: \n')
        fprintf('f:\n')
        disp(f_s)
        fprintf('h:\n')
        disp(h_s)
    end
    %% REDEFINITION OF VARIABLES, SAVE THE MODEL AND START AGAIN (IF IT IS NECESSARY)
    p(cl-txw)=[];   %   Remove parameters of the vector
    f=f_s;
    h=h_s;
    % Save the new model
    resultsname = sprintf('New_Model',modelname,date);
    fullresultsname = strcat(nmf,filesep,'results',filesep,resultsname);
    save(fullresultsname,'x','p','u','w','h','f','ics','known_ics');
    if it<num_repar
        
        %   There are still pending repairs
        STRIKE_GOLDD('New_Model',1);
        modelname='New_Model';
        [~,paths,opts,submodels,prev_ident_pars] = options();
        resultsname = sprintf('id_results_%s_%s',modelname,date);
        load(resultsname);
        fprintf(' -------------------------------------------------- \n');
        fprintf('Searching for symmetries of %s...\n', modelname);
        [transf,nuevas_variables,allVar]=Lie_Symmetry('New_Model');      
        
        % Verify that the transformation vector is not empty
        if (isempty(transf)==1)
           return 
        end
    else
        %   Reparametrized model
        fprintf('\n------------------------------------------------------------------\n');
        fprintf('>>> The model reparameterization has been completed successfully \n');
        fprintf('------------------------------------------------------------------\n');
        fprintf('The new FISPO model is: \n')
        fprintf('f:\n')
        disp(f)
        fprintf('h:\n')
        disp(h)
    end
end



end