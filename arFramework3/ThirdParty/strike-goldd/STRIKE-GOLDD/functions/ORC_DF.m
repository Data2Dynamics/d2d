%=========================================================================%
%====ORC-DF: An implementation of the Observability Rank Criterion for====%
%=====================systems with Direct Feedthrough=====================%

function ORC_DF(modelname,opts,prev_ident_pars)

tStart=tic;

%============================Load model===================================%

%(Optional) Create multi-experiment model:
if opts.multiexp
    ME_analysis(modelname,opts);
    modelname=strcat(modelname,'_',num2str(opts.multiexp_numexp),'Exp');
end

affine_modelname = sprintf('affine_%s',modelname);
affine_model_file=strcat(pwd,filesep,'models',filesep,strcat(affine_modelname,'.mat'));
exist_affine_model=exist(affine_model_file,'file');

%If the model has already been converted into affine in inputs form:
if exist_affine_model==2
    modelname=affine_modelname;
end

load(modelname) %#ok<*LOAD>

fprintf('\n >>> Analyzing observability of %s with affine in control system algorithm\n',modelname);

%======================Dimensions of the problem==========================%

%Number of states:
ns=numel(x);   %#ok<*USENS>
%Number of outputs:
m=numel(h);
%Number of unknown parameters:
if exist('p','var')
    np=numel(p);
else
    p=[];
    np=0;
end
%Number of unknown inputs:
if exist('w','var')&&numel(w)>0 %#ok<*NODEF>
    nw=numel(w);
    if opts.multiexp
        if opts.multiexp_user_nnzDerW
            opts.nnzDerW=opts.multiexp_nnzDerW;
        else
            opts.nnzDerW=repmat(opts.nnzDerW,1,opts.multiexp_numexp);
        end
    end
    if exist_affine_model~=2
    %Chek that system is affine in unknown inputs:    
        syms coeff_w coefh_w
        for j=1:nw
            for i=1:ns
                try
                    coeff_w=coeffs(f(i),w(j),'all');
                    if numel(coeff_w)>2
                        warning('Unable to convert model to affine in inputs form. System dynamics is not affine in unknown inputs.')
                        return
                    end
                catch
                    warning('Unable to convert model to affine in inputs form. System dynamics is not affine in unknown inputs.')
                    return
                end
                variables=symvar(coeff_w);
                nz_variables=simplify(variables-w(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. System dynamics is not affine in unknown inputs.')
                    return
                end
            end
            for i=1:m
                try 
                    coefh_w=coeffs(h(i),w(j),'all');
                    if numel(coefh_w)>2
                        warning('Unable to convert model to affine in inputs form. Output dynamics is not affine in unknown inputs.')
                        return
                    end
                catch
                    warning('Unable to convert model to affine in inputs form. Output dynamics is not affine in unknown inputs.')
                    return
                end
                variables=symvar(coefh_w);
                nz_variables=simplify(variables-w(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in unknown inputs.')
                    return
                end
            end
        end
    end

    %Create array of unknown inputs derivatives and set certain derivatives to zero:
    w_der=sym(zeros(nw,opts.affine_kmax+1));
    if numel(opts.nnzDerW)==nw
      for ind_w=1:nw
          if opts.nnzDerW(ind_w)<opts.affine_kmax+1
              w_der(ind_w,1:opts.nnzDerW(ind_w))=sym(strcat(char(w(ind_w)),sprintf('d')),[1 opts.nnzDerW(ind_w)]);
          else
              w_der(ind_w,:)=sym(strcat(char(w(ind_w)),sprintf('d')),[1 opts.affine_kmax+1]);
          end
      end   
    elseif numel(opts.nnzDerW)>0
      warning('The size of vector that contains the number of nonzero unknown inputs derivatives must be %d',nw)
      return
    else
    %if opts.nnzDerw=[],it is assumed that opts.nnzDerw(i)=inf, i=1,...,nw:
      for ind_w=1:nw, w_der(ind_w,:)=sym(strcat(char(w(ind_w)),sprintf('d_')),[1 opts.affine_kmax+1]);end
    end   
else
    w=[];
    nw=0;
    w_der=sym(zeros(nw,opts.affine_kmax+1));
end

%Number of known inputs:
if exist('u','var')&& numel(u)>0
    nu=numel(u);
    %Convert system into affine in inputs form:
    if exist_affine_model~=2
        fprintf('\n >>> Converting control system into affine in inputs form...\n');
        syms coeff_u coefh_u
        %Initialize known inputs coefficients of system dynamics: 
        fu=sym(zeros(ns,nu));
        %Initialize known inputs coefficients of outputs dynamics:
        hu=sym(zeros(m,nu));
        %Chek that system is affine in known inputs and almacenate coefficients:
        for j=1:nu
            for i=1:ns
                try
                    coeff_u=coeffs(f(i),u(j),'all');
                    if numel(coeff_u)>2
                        warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                        return
                    end
                catch
                    warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                    return
                end
                variables=symvar(coeff_u);
                nz_variables=simplify(variables-u(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. System dynamics is not affine in known inputs.')
                    return
                end
                %Almacenate uk coefficients of system dynamics: 
                if numel(coeff_u)==2, fu(i,j)=coeff_u(1);end
            end
            for i=1:m
                try
                    coefh_u=coeffs(h(i),u(j),'all');
                    if numel(coefh_u)>2
                        warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                        return
                    end
                catch
                    warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                    return
                end
                variables=symvar(coefh_u);
                nz_variables=simplify(variables-u(j)*ones(1,numel(variables)));
                if numel(find(nz_variables==0))~=0
                    warning('Unable to convert system to affine in inputs form. Output dynamics is not affine in known inputs.')
                    return
                end
                %Almacenate uk coefficients of output dynamics:
                if numel(coefh_u)==2, hu(i,j)=coefh_u(1);end
            end
        end
        %State-unknown inputs 0-augmented system dynamics:
        fxw=simplify(f-fu*u);
        %State-unknown inputs output dynamics:
        hxw=simplify(h-hu*u);
        %Remove dependent columns of fu:
        rank_fu=rank(fu);
        if rank_fu<min(ns,nu)
            zc_fu=[];
            for j=1:nu
                aux_fu=fu;
                aux_fu(:,j)=[];
                if rank(aux_fu)==rank_fu
                    zc_fu=[zc_fu j];
                end
            end
            fu(:,zc_fu)=[];
        end
        %Remove dependent columns of hu:
        rank_hu=rank(hu);
        if rank_hu<min(m,nu)
            zc_hu=[];
            for j=1:nu
                aux_hu=fu;
                aux_hu(:,j)=[];
                if rank(aux_hu)==rank_hu
                    zc_hu=[zc_hu j];
                end
            end
            hu(:,zc_hu)=[];
        end
        % Save affine model:
        fullaffinename = strcat(pwd,filesep,'models',filesep,affine_modelname);
        if exist('ics','var')==0,ics=[];end
        if exist('known_ics','var')==0,known_ics=[];end
        save(fullaffinename,'x','u','w','p','f','h','fu','fxw','hu','hxw','ics','known_ics');
    end
else
    u=[];
    nu=0;
    fu=[];
    hu=[];
    fxw=f;
    hxw=h;
end

%Remove parameters that have already been classified as identificable:
if exist('prev_ident_pars','var')
    for i=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(i),p);
        p(original_index)=[]; %#ok<*AGROW>
    end
    np=numel(p);
end

fprintf('\n >>> The model contains:\n %d states:\n %s',ns,char(x));
fprintf('\n %d outputs:\n %s',m,char(h));
fprintf('\n %d known inputs:\n %s',nu,char(u));
fprintf('\n %d unknown inputs:\n %s',nw,char(w));
fprintf('\n %d parameters:\n %s',np,char(p));

%====================(Optional)Parallel preferences=======================%

if opts.affine_parallel
    fprintf('\n >>> Initializating parallel preferences...')
    pool=gcp('nocreate');
    mycluster=parcluster('local');
    max_workers=mycluster.NumWorkers;
    if isempty(pool)
        if opts.affine_workers<=max_workers
            parpool(opts.affine_workers);
        else
            warning('The maximum number of workers is %d. Introduce a number of workers up to %d.',max_workers,max_workers)
            return
        end
    elseif pool.NumWorkers~=opts.affine_workers
        if opts.affine_workers<=max_workers
            delete(pool);
            parpool(opts.affine_workers);
        else
            warning('The maximum number of workers is %d. Introduce a number of workers up to %d.',max_workers,max_workers)
            return
        end
    end
end

%%===========================Initialization==============================%%

%Initialize computation time of each stage:
stage_time=zeros(opts.affine_kmax+1,1);
%Stage_time can be separated (essentially) in:
%Time spent calculating the differential of Omega:
dif_time=zeros(opts.affine_kmax+1,1);
%Time spent calculating the rank of Omega:
rank_time=zeros(opts.affine_kmax+1,1);
%Time spent calculating partial ranks of Omega:
partial_rank_time=zeros(opts.affine_kmax+1,1);

%Initialize number of states in each stage:
nxau=zeros(1,opts.affine_kmax+1);
%Number of states of 0-augmented system:
nxau(1)=ns+np+nw;
%Maximum number of states:
max_ns=nxau(1)+numel(find(w_der));
%Initialize rank of Omega for each stage:
rango=zeros(1,opts.affine_kmax);
%Initialize matrix whose components are the values of partial ranks for each state and iteration:
partial_ranks=zeros(max_ns,opts.affine_kmax);

tStage=tic;
%Stage counter:
k=0;
%0-augmented system state:
xaug=[x;p;w];
%Reshape hu as a column vector:
hu=reshape(hu,[],1); 
%System dynamics of 0-augmented system:
fxw=[fxw;zeros(np+nw,1)]+[zeros(ns+np,1);w_der(:,1)]; 
%Known inputs contribution of 0-augmented system dynamics:
if numel(fu)~=0
    fu=[fu;zeros(np+nw,numel(fu(1,:)))]; 
end
%Initialize matrix which (i,k)-entry is =0 if the ith state is obs. at the kth stage and =1 if it is not:
unobs_states_ind=ones(max_ns,opts.affine_kmax);
%First codistribution vector:
Delta_Omega=[hxw;hu];
tDer=tic;
%First derivative of codistribution vector:
dif_Delta_omega=jacobian(Delta_Omega,xaug);
dif_time(1)=toc(tDer);
%Initialize observability matrix:
dif_Omega=dif_Delta_omega;
%Index actualization:
k=k+1;
%0-stage computation time:
stage_time(1)=toc(tStage);

%(Optional) Define numerical equivalents of the symbolic variables:
if opts.numeric == 1
    [~,~,nz_w_der]=find(w_der);
    xau=[xaug;[f;zeros(np,1);reshape(nz_w_der,[],1)];h;u];
    allvariables = symvar(xau);
    numbers = vpa(0.1+rand(size(allvariables)));
end  

%==============================Kth stage==================================%

stage_time(2)=stage_time(1);

while k<=opts.affine_kmax && stage_time(k+1)-stage_time(k)<opts.affine_tStage

    fprintf('\n >>> Building observability matrix of %d-augmented system...',k)

    tStage=tic;

    %Actualization of delta_omega (Lie derivative of Omega): 
    Delta_Omega=dif_Delta_omega*[fxw fu]; 
    %Reshape delta_omega as a column vector:
    Delta_Omega=reshape(Delta_Omega,[],1);
    %Including as states only nonzero kth derivatives of unknown inputs:
    [nz_r,~,nz_w_der] = find(w_der(:,k)); 
    xaug=[xaug;nz_w_der];
    %Number of states of k-augmented system:
    nxau(k+1)=numel(xaug);
    %Calculate the differential of Omega:
    tDer=tic;
    %Calculate the differential of Delta_Omega:
    dif_Delta_omega=jacobian(Delta_Omega,xaug);
    %Actualization of the differential of Omega:
    dif_Omega=[dif_Omega zeros(numel(dif_Omega(:,1)),nxau(k+1)-nxau(k)); dif_Delta_omega];
    %Computation time of the differential of Omega:
    dif_time(k+1)=dif_time(k)+toc(tDer);

%==============Investigate observability of k-augmented system============%

    fprintf('\n >>> Calculating rank of %d-augmented system observability matrix...',k)

    %(Optional) Replace known initial conditions:
    if opts.replaceICs == 1
        xind = find(known_ics);
        if size(ics) ~= size(x)
            ics = transpose(ics);
        end
        dif_Omega= subs(dif_Omega,x(xind),ics(xind));
    end 

    %(Optional) Replace numerical equivalents of symbolic variables:
    if opts.numeric == 1
        num_dif_omega = subs(dif_Omega,allvariables,numbers); 
    else
        num_dif_omega = dif_Omega;
    end 

    %Calculate rank of k-observability matrix:
    tRank=tic;
    rango(k)=rank(num_dif_omega);
    %Actualization of rank computation time:
    last_rank_time=toc(tRank);
    rank_time(k)=rank_time(k)+last_rank_time;
    rank_time(k+1)=rank_time(k);

    fprintf('\n     Rank = %d (calculated in %d seconds)',rango(k),last_rank_time)

    %Update partial ranks matrix:
    partial_ranks(1:nxau(k+1),k)=rango(k)-1;
    %Store indices of unobservable states:
    [nf_unobs,~]=find(unobs_states_ind(1:nxau(k+1),k));

    %Check if the model is FISPO:
    if rango(k)==nxau(k+1)
        %Print results if the model is FISPO:
        stage_time(k+1)=stage_time(k)+toc(tStage);

        fprintf('\n\n ------------------------ \n');
        fprintf(' >>> RESULTS SUMMARY:\n');
        fprintf(' ------------------------ \n');

        fprintf('\n >>> The model is k-row observable for k = %d  \n',k)
        fprintf('\n >>> The model is Fully Input-State-Parameter Observable (FISPO):');
        if nw>0, fprintf('\n All its unknown inputs are observable.'); end
        fprintf('\n All its states are locally structurally observable.');
        if np>0, fprintf('\n All its parameters are locally structurally identifiable.'); end

        if opts.affine_graphics
        %All states are k-row observable:
        unobs_states_ind(:,k)=zeros(max_ns,1);
        %Ticks including system states for labelling y-axis:
        xau_ticks=flip(arrayfun(@char, xaug, 'uniform',0));                                     
        for i=1:k
            %Unobservable states at k=i:
            unobs_i=find(unobs_states_ind(1:nxau(i+1),i));
            %Number of unobservable states at k=i:
            n_unobs=numel(unobs_i);
            %Plot unobservable and observable states for {1,...,k}:            
            figure(1)
            %Print x-tick at unobservable states:
            scatter(i*ones(n_unobs,1),(nxau(k+1)+1)*ones(n_unobs,1)-unobs_i,50,'x','b','LineWidth',1)
            grid on
            hold on
            %Observable states at k=i:
            obs_i=find(unobs_states_ind(1:nxau(i+1),i)==0);  
            %Number of observable states at k=i:
            n_obs=numel(obs_i);
            %Print o-tick at observable states:
            scatter(i*ones(n_obs,1),(nxau(k+1)+1)*ones(n_obs,1)-obs_i,50,'o','b','LineWidth',1)
            %Undefined states at k=i:
            non_def_states=(nxau(i+1)+1):nxau(k+1);
            %Print black point at undefined states:
            scatter(i*ones(1,nxau(k+1)-nxau(i+1)),(nxau(k+1)+1)*ones(1,nxau(k+1)-nxau(i+1))-non_def_states,70,'.','k','LineWidth',1.5)
        end
        legend('Unobservable states','Observable states','Non-defined states','location','southoutside','Orientation','horizontal');
        axis([0 k+1 0 nxau(k+1)+1]);                        
        xticks(1:k);
        xlabel('Stage');
        yticks(1:nxau(k+1));
        yticklabels(xau_ticks);
        title('Classification of model variables')
        hold off

        %Plot computation times for each iteration:
        dif_time=dif_time(1:k+1);
        stage_time=stage_time(1:k+1);
        rank_time=rank_time(1:k);
        partial_rank_time=partial_rank_time(1:k);

        figure(2)
        subplot(1,4,1)
        semilogy(0:k,dif_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xlabel('Stage')
        ylabel('t[s]')
        xticks(0:k)
        yticks(dif_time)
        ytickformat('%.3f')
        xlim([0 k])
        ylim([dif_time(1) dif_time(end)])
        title('Partial derivatives computation time')

        subplot(1,4,2)
        semilogy(1:k,rank_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xlabel('Stage')
        ylabel('t[s]')
        xticks(1:k)
        yticks(rank_time)
        ytickformat('%.3f')
        xlim([1 k])
        ylim([rank_time(1) rank_time(end)])
        title('Rank computation time')

        subplot(1,4,3)
        semilogy(1:k,partial_rank_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xlabel('Stage')
        ylabel('t[s]')
        xticks(1:k)
        yticks(partial_rank_time(1:k-1))
        ytickformat('%.3f')
        xlim([1 k])
        ylim([partial_rank_time(1) partial_rank_time(end)])
        title('Partial ranks computation time')

        subplot(1,4,4)
        semilogy(0:k,stage_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xlabel('Stage')
        ylabel('t[s]')
        xticks(0:k)
        yticks(stage_time)
        ytickformat('%.3f')
        xlim([0 k])
        ylim([stage_time(1) stage_time(end)])
        title('Stage computation time')
        end

        totaltime = toc(tStart);
        fprintf('\n Total execution time: %d \n\n',totaltime);

        %Save results if model is FISPO:
        resultsname = sprintf('id_results_%s',modelname);
        fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname,'_',date);
        warning off 'parallel:cluster:CannotSaveCorrectly'
        save(fullresultsname);

        % Delete affine model if necessary:
        if opts.affine_delete_model, delete(affine_model_file); end

        return 

    elseif nw==0 && k>1 && rango(k)==rango(k-1)

        %If rank(Omega_k)=rank(Omega_k-1), the rank will not increase anymore and system is not observable:
        stage_time(k+1)=stage_time(k)+toc(tStage);

        fprintf('\n >>> Observability matrix dimension will not increase including further derivatives.')
        fprintf('\n >>> The control system is not k-row observable for k higher or equal to 1.')

        if opts.affine_graphics
        %Ticks including system states for labelling y-axis:
        xau_ticks=flip(arrayfun(@char, xaug, 'uniform',0)); 
        for i=1:k
            %Unobservable states at k=i:
            unobs_i=find(unobs_states_ind(1:nxau(i+1),i));
            %Number of unobservable states at k=i:
            n_unobs=numel(unobs_i);
            %Plot unobservable and observable states for {1,...,k}:            
            figure(1)
            %Print x-tick at unobservable states:
            scatter(i*ones(n_unobs,1),(nxau(k+1)+1)*ones(n_unobs,1)-unobs_i,50,'x','b','LineWidth',1)
            %Plot settings:
            axis([0 k+1 0 nxau(k+1)+1]);                        
            xticks(1:k);
            xlabel('Stage');
            yticks(1:nxau(k+1));
            yticklabels(xau_ticks);
            grid on
            hold on
            %Observable states at k=i:
            obs_i=find(unobs_states_ind(1:nxau(i+1),i)==0);  
            %Number of observable states at k=i:
            n_obs=numel(obs_i);
            %Print o-tick at observable states:
            scatter(i*ones(n_obs,1),(nxau(k+1)+1)*ones(n_obs,1)-obs_i,50,'o','b','LineWidth',1)
            %Undefined states at k=i:
            non_def_states=(nxau(i+1)+1):nxau(k+1);
            %Print black point at undefined states:
            scatter(i*ones(1,nxau(k+1)-nxau(i+1)),(nxau(k+1)+1)*ones(1,nxau(k+1)-nxau(i+1))-non_def_states,70,'.','k','LineWidth',1.5)
        end
        legend('Unobservable states.','Observable states.','Non-defined states.','location','southoutside','Orientation','horizontal')
        hold off
        title('Classification of model variables')

        dif_time=dif_time(1:k+1);
        stage_time=stage_time(1:k+1);
        rank_time=rank_time(1:k);
        partial_rank_time=partial_rank_time(1:k);

        figure(2)
        subplot(1,4,1)
        semilogy(0:k,dif_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xlim([0 k])
        ylim([dif_time(1) dif_time(end)])
        xticks(0:k)
        yticks(dif_time)
        ytickformat('%.3f')
        xlabel('Stage')
        ylabel('t[s]')
        title('Partial derivatives computation time')

        subplot(1,4,2)
        semilogy(1:k,rank_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xlim([1 k])
        ylim([rank_time(1) rank_time(end)])
        xticks(1:k)
        yticks(rank_time)
        ytickformat('%.3f')
        xlabel('Stage')
        ylabel('t[s]')
        title('Rank computation time')

        subplot(1,4,3)
        semilogy(1:k,partial_rank_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
        grid on
        xticks(1:k)
        yticks(partial_rank_time(1:k-1))
        ytickformat('%.3f')
        xlabel('Stage')
        ylabel('t[s]')
        xlim([1 k])
        ylim([partial_rank_time(1) partial_rank_time(end)])
        title('Partial ranks computation time')

        subplot(1,4,4)
        semilogy(0:k,stage_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',10)
        grid on
        xlim([0 k])
        ylim([stage_time(1) stage_time(end)])
        xticks(0:k)
        yticks(stage_time)
        ytickformat('%.3f')
        xlabel('Stage')
        ylabel('t[s]')
        title('Stage computation time')
        end    

        %Observable and unobservable states:
        [nf_obs_states,~]=find(unobs_states_ind(1:ns,k)==0);
        obs_states=xaug(nf_obs_states);

        [nf_unobs_states,~]=find(unobs_states_ind(1:ns,k));
        unobs_states=xaug(nf_unobs_states);

        unobs_states_ind(1:ns,k)=1;

        %Identifiable and unidentifiable parameters:
        [nf_obs_par,~]=find(unobs_states_ind(1:(ns+np),k)==0);
        obs_par=xaug(nf_obs_par);

        unobs_states_ind(1:ns,k)=0;

        [nf_unobs_par,~]=find(unobs_states_ind(1:(ns+np),k));
        unobs_par=xaug(nf_unobs_par);

        fprintf('\n\n ------------------------ \n');
        fprintf(' >>> RESULTS SUMMARY:\n');
        fprintf(' ------------------------ \n');

        if numel(obs_states)>0, fprintf('\n >>> The original observable states are: \n   %s', char(obs_states));end
        if numel(unobs_states)>0, fprintf('\n >>> The original unobservable states are: \n   %s', char(unobs_states));end
        if numel(obs_par)>0, fprintf('\n >>> The identifiable parameters are: \n    %s',char(obs_par)); end
        if numel(unobs_par)>0, fprintf('\n >>> The unidentifiable parameters are: \n    %s',char(unobs_par)); end

        totaltime=toc(tStart);
        fprintf('\n Total execution time: %d \n\n',totaltime);

        %Save results:
        resultsname = sprintf('id_results_%s',modelname);
        fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname,'_',date);
        warning off 'parallel:cluster:CannotSaveCorrectly'
        save(fullresultsname);

        % Delete affine model if necessary:
        if opts.affine_delete_model, delete(affine_model_file);end
        return  
    end

    fprintf('\n >>> The %d-augmented system is not FISPO.',k)

%=========Investigate partial observability of k-augmented system=========%

    fprintf('\n >>> Investigating partial observability of %d -augmented system...',k)

    tPartial=tic;
    if opts.affine_parallel
        %Partial rank computation (using parallel toolbox):
        %For the (k-1)-row unobs.states, remove its column and recalculate rank: 
        parfor i=1:numel(nf_unobs)                  
            elim_dif_omega=num_dif_omega;
            elim_dif_omega(:,nf_unobs(i))=[];     
            new_rank(i)=rank(elim_dif_omega);
        end
        %Update obs-unobs states array:
        for i=1:numel(nf_unobs)
            if new_rank(i)<rango(k)
                unobs_states_ind(nf_unobs(i),k)=0;
                partial_ranks(nf_unobs(i),k)=new_rank(i);
            end
        end
        partial_rank_time(k)=partial_rank_time(k)+toc(tPartial);
        partial_rank_time(k+1)=partial_rank_time(k);    
    else
        %Partial rank calculation (without parallel toolbox):
        for i=1:numel(nf_unobs)
            elim_dif_omega=num_dif_omega; 
            elim_dif_omega(:,nf_unobs(i))=[];      
            partial_ranks(nf_unobs(i),k)=rank(elim_dif_omega);
            if partial_ranks(nf_unobs(i),k)<rango(k), unobs_states_ind(nf_unobs(i),k)=0;end
        end
        partial_rank_time(k)=partial_rank_time(k)+toc(tPartial);
        partial_rank_time(k+1)=partial_rank_time(k);
    end

%==================Actualization of system dynamics=======================%

    %Store unobs-obs. indices at the k-th stage:
    unobs_states_ind(:,k+1)=unobs_states_ind(:,k); 
    %Update state-unknown inputs system dynamics:
    fxw=[fxw;w_der(nz_r,k+1)];
    %Update known inputs contribution to system dynamics:
    fu=[fu;zeros(nxau(k+1)-nxau(k),nu)];
    %Update stage computation time:
    stage_time(k+1)=stage_time(k)+toc(tStage);
    %Update index:
    k=k+1;
end 

%================Observability results if k=opts.kmax+1===================%

if opts.affine_graphics
    %Ticks with states for labelling y-axis:
    xau_ticks=flip(arrayfun(@char, xaug, 'uniform',0));
    %Plot unobservable and observable states for k=1,...,opts.kmax:
    figure(1)
    for i=1:k-1
        %Unobservable states at the ith stage:
        unobs_i=find(unobs_states_ind(1:nxau(i+1),i));     
        %Number of unobservable states at the ith stage:
        n_unobs=numel(unobs_i);
        %Print x-tick at unobservable states:
        scatter(i*ones(n_unobs,1),(nxau(k)+1)*ones(n_unobs,1)-unobs_i,50,'x','b','LineWidth',1)
        grid on
        hold on
        %Observable states at the ith stage:
        obs_i=find(unobs_states_ind(1:nxau(i+1),i)==0);  
        %Number of observable states at the ith stage:
        n_obs=numel(obs_i);
        %Print o-tick at observable states:
        scatter(i*ones(n_obs,1),(nxau(k)+1)*ones(n_obs,1)-obs_i,50,'o','b','LineWidth',1)
        %Undefined states at the ith stage:
        non_def_states=(nxau(i+1)+1):nxau(k);
        %Print black dot at undefined states:
        scatter(i*ones(1,nxau(k)-nxau(i+1)),(nxau(k)+1)*ones(1,nxau(k)-nxau(i+1))-non_def_states,70,'.','k','LineWidth',1.5)
    end
    legend('Unobservable states.','Observable states.','Non-defined states.','Location','southoutside','Orientation','horizontal')
    hold off
    title('Classification of model variables')
    xlabel('Stage')
    axis([0 k 0 nxau(k)+1]);                        
    xticks(1:k-1);
    yticks(1:nxau(k));
    yticklabels(xau_ticks);

    %Plot computation times:

    dif_time=dif_time(1:k);
    stage_time=stage_time(1:k);
    rank_time=rank_time(1:k-1);
    partial_rank_time=partial_rank_time(1:k-1);

    figure(2)
    subplot(1,4,1)
    semilogy(0:k-1,dif_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
    grid on
    xlim([0 k-1])
    ylim([dif_time(1) dif_time(end)])
    xlabel('Stage')
    ylabel('t[s]')
    xticks(0:k-1)
    yticks(dif_time)
    ytickformat('%.3f')
    title('Partial derivatives computation time')

    subplot(1,4,2)
    semilogy(1:k-1,rank_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
    grid on
    xlabel('Stage')
    ylabel('t[s]')
    xticks(1:k-1)
    yticks(rank_time)
    xlim([1 k-1])
    ylim([rank_time(1) rank_time(end)])
    ytickformat('%.3f')
    title('Rank computation time')

    subplot(1,4,3)
    semilogy(1:k-1,partial_rank_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
    grid on
    xlabel('Stage')
    ylabel('t[s]')
    xticks(1:k-1)
    yticks(partial_rank_time)
    ytickformat('%.3f')
    xlim([1 k-1])
    ylim([partial_rank_time(1) partial_rank_time(end)])
    title('Partial ranks computation time')

    subplot(1,4,4)
    semilogy(0:k-1,stage_time,'LineStyle','--','Linewidth',1.5,'Color',[0.5 0.7 0.3],'Marker','^','MarkerSize',8)
    grid on
    xlabel('Stage')
    ylabel('t[s]')
    xticks(0:k-1)
    yticks(stage_time)
    xlim([0 k-1])
    ylim([stage_time(1) stage_time(end)])
    ytickformat('%.3f')
    title('Stage computation time')
end

if k>opts.affine_kmax
    fprintf('\n\n >>> Maximum number of iterations or rebased.');
    fprintf('\n >>> You can increase it by changing <<opts.affine_kmax>> (currently opts.affine_kmax = %d)',opts.affine_kmax)
else
    fprintf('\n\n >>> Maximum computation time of the last iteration rebased.');
    fprintf('\n >>> You can increase it by changing <<opts.affine_tStage>> (currently opts.affine_tStage = %d)',opts.affine_tStage)
end

fprintf('\n\n ------------------------ \n');
fprintf(' >>> RESULTS SUMMARY:\n');
fprintf(' ------------------------ \n');

% Store indices of observable states:
[nf_obs_states,~]=find(unobs_states_ind(1:ns,k)==0);
obs_states=xaug(nf_obs_states);

unobs_states_ind(1:ns,k)=1;

% Store indices of identifiable parameters:
[nf_obs_par,~]=find(unobs_states_ind(1:(ns+np),k)==0);
obs_par=xaug(nf_obs_par);

unobs_states_ind((ns+1):(ns+np),k)=1;

% Store indices of observable inputs:
if nw>0
    [nf_obs_inputs,~]=find(unobs_states_ind(:,k)==0);
    obs_inputs=xaug(nf_obs_inputs);
end

if numel(obs_states)>0, fprintf('\n >>> The original observable states are: \n   %s', char(obs_states));end
if numel(obs_par)>0, fprintf('\n >>> The observable parameters are: \n    %s',char(obs_par)); end
if nw>0 && numel(obs_inputs)>0, fprintf('\n >>> The observable unknown inputs are: \n   %s',char(obs_inputs));end
fprintf('\n\n >>> The model is not k-row observable for k < %d  \n',k);

totaltime=toc(tStart);
fprintf('\n\n Total execution time: %d \n\n',totaltime);

% Save results:
resultsname = sprintf('id_results_%s',modelname);
fullresultsname = strcat(pwd,filesep,'results',filesep,resultsname,'_',date);
warning off 'parallel:cluster:CannotSaveCorrectly'
save(fullresultsname);

% Delete affine model (optional):
if opts.affine_delete_model, delete(affine_model_file);end

end
