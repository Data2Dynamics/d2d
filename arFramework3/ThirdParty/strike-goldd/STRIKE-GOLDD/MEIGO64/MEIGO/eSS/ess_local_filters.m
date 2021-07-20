% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ess_local_filters.m 1133 2013-12-02 12:54:12Z attila $
if stage_1==1
    %In this first stage we apply the local search over the best solution
    %found so far
    [minimo,I]=min(child_values_penalty);
    if minimo<fbest
        x0=child(I,:);
        f0=minimo;
    else
        x0=xbest;
        f0=fbest;
    end

    %Switch to stage 2
    n_critico=local_n2;

    if local_bestx
        use_bestx=0;    %not to go directly to stage 2
    end
    %use_bestx is a flag indicating that the best solution has been
    %improved
elseif stage_2==1;
    if local_bestx && use_bestx
        x0=xbest;
        f0=fbest;
        use_bestx=0;
    else
        %Calculate distance between our children and local_solutions+initial points
        local_init=[local_solutions; initial_points];
        xuxl=x_U-x_L;
        child_norm=child./repmat(xuxl,size(child,1),1);
        local_init_norm=local_init./repmat(xuxl,size(local_init,1),1);
        dist=eucl_dist(child_norm',local_init_norm');
        distance=dist';
        [eee ooo]=sort(-min(distance,[],1));
        [aaa uuu]=sort(child_values_penalty);
        child_points=zeros(1,length(child_values_penalty));

        www=local_balance;
        for i=1:length(child_values_penalty)
            child_points(ooo(i))=child_points(ooo(i))+www*i;
            child_points(uuu(i))=child_points(uuu(i))+(1-www)*i;
        end
        [minimum,I]=min(child_points);
        x0=child(I,:);
        f0=child_values_penalty(I);
    end
end

if iterprint
    fprintf('Call local solver: %s \n', upper(local_solver))
    fprintf('Initial point function value: %f \n',f0);
    tic
end

%%% AFV 19/04/17: try-catch to avoid crashes due to errors in local search:
try
    if isFjacDefined
            [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,fjac,...
                local_solver,local_iterprint,local_tol,weight,nconst,tolc,opts.local,varargin{:});
        else
            [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,[],...
                local_solver,local_iterprint,local_tol,weight,nconst,tolc,opts.local,varargin{:});
    end
catch
    fprintf('The local search crashed \n');
    x        = x0;
    fval     = NaN;
    exitflag = -1;
    numeval  = 1;
end
%%% AFV: end of changes 19/04/17

if iterprint
    fprintf('Local solution function value: %g \n',fval);
    fprintf('Number of function evaluations in the local search: %i \n',numeval);
    fprintf('CPU Time of the local search: %f  seconds \n\n',toc);
end

initial_points=[initial_points; x0];

%Careful with solnp and the reported number of evaluations

nfuneval=nfuneval+numeval;
if stage_1==1
    stage_1=0;
    stage_2=1;
end

%Evaluate the local solution
[val,val_penalty,pena,nlc,include,x] = ssm_evalfc(x,x_L,x_U,fobj,nconst,c_L,c_U,...
    tolc,weight,int_var,bin_var,varargin{:});
nfuneval=nfuneval+1;
if include && pena<=tolc
    if val_penalty<fbest
        fbest=val_penalty;
        xbest=x;
        if fbest<fbest_lastiter
            Results.f=[Results.f fbest];
            Results.x=[Results.x; xbest];
            Results.time=[Results.time toc(tinit)];
            Results.neval=[Results.neval nfuneval];
            if plot_results==1
                stairs(Results.time,Results.f)
                xlabel('CPU Time (s)');
                ylabel('Objective Function Value (-)');
                title(strcat('Best f:  ',num2str(fbest)))
                drawnow
            end
        end
    end
end

%Check the stopping criterion
if nfuneval>= maxeval
    fin=1;
    blacklist=ones(dim_refset);
elseif toc(tinit)>=maxtime
    fin=2;
    blacklist=ones(dim_refset);
elseif not(isempty(vtr))
    if fbest<=vtr
        fin=3;
        blacklist=ones(dim_refset);
    end
end

if include
    if not(isempty(local_solutions))
        %Add the found local solution to the list (if it is a new one)
        adicionar_local=1;
        [f,ind]=ssm_isdif2(x,local_solutions,1e-2,1);
        %If it is too close to a previously found local solution, we do not
        %add it
        if f
            adicionar_local=0;
        end
    else
        %Add it if no local solutions have been found yet
        adicionar_local=1;
    end
    if adicionar_local==1

        %choose either this (replace the worst element in Refset)

        %         [aaa,jjj]=max(Refset_values_penalty);
        %
        %         if val_penalty<Refset_values_penalty(jjj)
        %             Refset(jjj,:)=x;
        %             Refset_values(jjj)=val;
        %             Refset_values_penalty(jjj)=val_penalty;
        %             Refset_nlc(jjj,:)=nlc;
        %             penalty(jjj)=pena;
        %         end

        %or this (Replace the initial point)

        %         if val_penalty<f0
        %             child(I,:)=x;
        %             child_values(I)=val;
        %             child_values_penalty(I)=val_penalty;
        %             child_nlc(I,:)=nlc;
        %         end

        %or this (Replace the parent of the initial point)

        if val_penalty<Refset_values_penalty(child_parent_index(I))
            Refset(child_parent_index(I),:)=x;
            Refset_values(child_parent_index(I))=val;
            Refset_values_penalty(child_parent_index(I))=val_penalty;
            Refset_nlc(child_parent_index(I),:)=nlc;
            penalty(child_parent_index(I))=pena;

            refset_change(child_parent_index(I))=0;
        end
        local_solutions=[local_solutions;x];
        local_solutions_values=[local_solutions_values;fval];
    end
end

