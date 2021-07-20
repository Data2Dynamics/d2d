% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_localsolver.m 2031 2015-08-24 11:49:26Z attila $
function [x,fval,exitflag,numeval]=ssm_localsolver(X0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,fjac,...
    local_solver,local_iterprint,local_tol,weight,nconst,tolc,local_opts,varargin)
global ccll ccuu n_upper n_lower n_fun_eval

global hjfun hjxl hjxu hjcl hjcu hjweight hjtolc
hjfun=fobj;
hjxl=x_L;
hjxu=x_U;
hjcl=c_L;
hjcu=c_U;
hjweight=weight;
hjtolc=tolc;

nvar=length(X0);

%Scaling factors (can be useful for some local solvers). This could be
%defined just once in the main ssm.m
scale.x=x_U-x_L;

f1=feval(fobj,X0,varargin{:});
n_fun_eval=1;
if nconst
    g_u=c_U;
    g_l=c_L;
    
    g_l(abs(g_l)==inf)=1;
    g_u(abs(g_u)==inf)=1;
    
    g_l=abs(g_l);
    g_u=abs(g_u);
    
    scale.g=max(g_l,g_u);
else
    scale.g=1;
end

if abs(f1)==inf || isnan(f1)
    scale.f=1;
else
    scale.f=f1;
end






switch local_solver
    case 'clssolve'
        %clssolve
        %no acepta restricciones lineales de igualdad
        Prob = clsAssign(fobj, [], [], x_L, x_U, 'test', X0, ...
            [], [], 0, [],[],[],...
            [],[],[],nlc_fun,[],[],c_L,c_U);
        %Prob.optParam.MaxIter=50;
        
        if local_iterprint, Prob.PriLevOpt = 2; end
        
        switch local_tol
            case 1
                Prob.optParam.eps_f=tolc*100;
                Prob.optParam.eps_x=tolc*100;
                Prob.optParam.eps_c=tolc*100;
            case 2
                Prob.optParam.eps_f=tolc;
                Prob.optParam.eps_x=tolc;
                Prob.optParam.eps_c=tolc;
            case 3
                Prob.optParam.eps_f=tolc/100;
                Prob.optParam.eps_x=tolc/100;
                Prob.optParam.eps_c=tolc/100;
        end
        
        
        Result=tomrun('clssolve',Prob);                        %Aplicamos el solver
        
        
        x=Result.x_k';
        fval=Result.r_k;
        numeval=Result.ResEv;
        %if Result.Iter<Prob.optParam.MaxIter
        exitflag=1;
        %else
        %   exitflag=0;
        %end
        
    case 'snopt'
        %snopt
        
        %no acepta restricciones lineales de igualdad
        
        
        Prob = simAssign(fobj, [], [], [], x_L, x_U, [], X0, ...
            [], [], [], [], [], c_L, c_U);
        
        
        
        %                             Prob = conAssign(fobj, [], [], [],x_L, x_U, 'Test', X0,...
        %         [],[],[],[],[],nlc_fun,[],[],[],c_L,c_U);
        %
        Prob.NumDiff=6;
        %Prob.optParam.MaxIter=250;
        
        if local_iterprint, Prob.PriLevOpt = 2; end;
        
        switch local_tol
            case 1
                Prob.optParam.eps_f=tolc*100;
                Prob.optParam.eps_x=tolc*100;
                Prob.optParam.eps_c=tolc*100;
            case 2
                Prob.optParam.eps_f=tolc;
                Prob.optParam.eps_x=tolc;
                Prob.optParam.eps_c=tolc;
            case 3
                Prob.optParam.eps_f=tolc/100;
                Prob.optParam.eps_x=tolc/100;
                Prob.optParam.eps_c=tolc/100;
        end
        
        Result=tomrun('snopt',Prob);                        %Aplicamos el solver
        
        x=Result.x_k';
        fval=Result.f_k;
        numeval=Result.FuncEv;
        
        %if Result.MinorIter<Prob.optParam.MaxIter
        exitflag=1;
        %else
        %exitflag=0;
        %end
        
        
    case 'nomad'
        %nomadm
        global params
        params=varargin;
        nvar=length(X0);
        problem_path=pwd;
        %Fichero con el punto inicial
        fid=fopen('fobj_nomad_x0.m','w');
        fprintf(fid,'function iterate= fobj_nomad_x0\n');
        fprintf(fid,'iterate(1).x= [');
        for i=1:size(X0,2)
            fprintf(fid,' %g;',X0(i));
        end
        fprintf(fid,' ];\n');
        fprintf(fid,'iterate(1).p= {};\n');
        fclose(fid);
        
        switch local_tol
            case 1
                tol_nomad=tolc*100;
            case {2,3}
                tol_nomad=tolc*10;
                %             case 3
                %                 tol_nomad=tolc;
        end
        if local_iterprint
            close;
        end
        [BestF,BestI,RunStats,RunSet]=mads_batch(problem_path,tol_nomad,nvar,length(n_lower)+length(n_upper),local_iterprint);
        clear('fobj_nomad_x0.m');
        
        
        if ~isempty(BestF)
            x=BestF.x';
            fval=BestF.f;
            exitflag=1;
        else
            x=X0;
            fval=inf;
            exitflag=-1;
        end
        
        
        numeval=RunStats.nFunc;
        %if numeval<1000
        
        %else
        %exitflag=0;
        %end
        %keyboard
        
        
    case 'npsol'
        %npsol
        %no acepta restricciones lineales de igualdad
        
        Prob = simAssign(fobj, [], [], [], x_L, x_U, [], X0, ...
            [], [], [], [], [], c_L, c_U);
        %     Prob = conAssign, [], [], [],x_L, x_U, 'Test', X0,...
        %         [],[],[],[],[],nlc_fun,[],[],[],c_L,c_U);
        
        %Prob.optParam.MaxIter=250;
        
        %Prob.PriLevOpt = 2;
        Prob.NumDiff=6;
        if local_iterprint, Prob.PriLevOpt = 2; end;
        
        switch local_tol
            case 1
                Prob.optParam.eps_f=tolc*100;
                Prob.optParam.eps_x=tolc*100;
                Prob.optParam.eps_c=tolc*100;
            case 2
                Prob.optParam.eps_f=tolc;
                Prob.optParam.eps_x=tolc;
                Prob.optParam.eps_c=tolc
            case 3
                Prob.optParam.eps_f=tolc/100;
                Prob.optParam.eps_x=tolc/100;
                Prob.optParam.eps_c=tolc/100;
        end
        
        
        Result=tomrun('npsol',Prob);                        %Aplicamos el solver
        
        x=Result.x_k';
        fval=Result.f_k;
        numeval=Result.FuncEv;
        
        %if Result.MinorIter<Prob.optParam.MaxIter
        exitflag=1;
        %else
        %exitflag=0;
        %end
        
        
        
    case 'solnp'
        %solnp
        %Fichero con el punto inicial
        
        cl=c_L;
        cu=c_U;
        
        cl(1:neq)=[];
        cu(1:neq)=[];
        
        aaa=find(X0==x_L);
        bbb=find(X0==x_U);
        X0(aaa)=X0(aaa)+0.001*abs(X0(aaa))+1e-10;
        X0(bbb)=X0(bbb)-0.001*abs(X0(bbb))-1e-10;
        
        xxb=[X0' x_L' x_U'];
        
        lh=cl';
        uh=cu';
        
        
        lh(find(lh==-inf))=-1e5;
        uh(find(uh==inf))=1e5;
        
        ib=[lh uh];
        
        if isempty(ib)
            ib=0;
        end
        ro=1;
        %Parece que aumentando las iteraciones menores afina mas, aunque
        %presumiblemente va mas lento. Iria mejor para las iteraciones finales
        maj_it=10*length(X0);
        min_it=10;
        
        if local_tol==1
            delta=1e-5;
            epsilon=tolc*10;
        elseif local_tol==2
            delta=1e-5;
            epsilon=tolc;
        elseif local_tol==3
            delta=1e-5;
            epsilon=tolc/10;
        end
        
        
        
        op=[ro maj_it min_it 0 delta epsilon];
        
        %         [x,oh]=solnp(xxb,ib,op);
        %         x=x';
        %         fval=oh(end);
        %         numeval=n_fun_eval;
        %         exitflag=1;
        %
        
        [y1,oh,flag,inc,lambda,hess]=solnp(@snpfun,xxb,ib,op,[],[],local_iterprint,fobj,scale,nconst,varargin{:});
        x=y1';
        fval=oh(end);
        numeval=n_fun_eval;
        exitflag=flag;
        
    case 'nl2sol'
        % nl2sol of OPTI Toolbox
        
        % gradient calculation of the residual:
        if isempty(fjac)
            % if no Jacobian is comming from AMIGO but
            % ess_options orders to use mkl
            if strcmpi(local_opts.nl2sol.grad,'mkl')
                fjac = @(x) mklJac(@objf_nl2sol,x,ndata);
                local_opts.nl2sol = rmfield(local_opts.nl2sol,'grad');
            else
                % neither the jacobian is given and the default is not mkl
                fjac = [];
                % NL2sol will use internal Finite Differences
            end
        else
            % we already wrote the file fjac_nl2sol containing the
            % interface to the external jacobian (AMIGO_PEjac)
            fjac = @fjac_nl2sol;
        end
        
        [x, fval,exitflag]=feval(@call_nl2sol,X0,fjac, ndata, x_L, x_U,local_opts.nl2sol);
        
        numeval=n_fun_eval+1;
        
    case {'n2fb','dn2fb'}
        %n2fb
        x=feval(@call_n2fb,X0, ndata, length(x_U), x_L, x_U,local_solver);
        fval=feval(fobj,x,varargin{:});
        exitflag=1;
        numeval=n_fun_eval+1;
        
    case 'lbfgsb'
        % lbfgsb of OPTI Toolbox
        
        % gradient calculation of the residual:
        if isempty(fjac)
            if strcmpi(local_opts.lbfgsb.grad,'mkl')
                
                fjac = @(x) mklJac(@objf_lbfgsb,x,1);
                local_opts.lbfgsb = rmfield(local_opts.lbfgsb,'grad');
            else
                fjac = [];
            end
        else
            fjac = @fjac_lbfgsb;
        end
        
        [x, fval,exitflag]=feval(@call_lbfgsb,X0,fjac, x_L, x_U,local_opts.lbfgsb);
        
        numeval=n_fun_eval+1;
    case 'dhc'
        %dhc
        nvar=numel(X0);
        X0=(X0-x_L)./(x_U-x_L);
        
        %initsize=max((x_U-x_L)/2);
        initsize=0.1;
        
        
        if local_tol==1
            thres=1e-6;
        elseif local_tol==2
            thres=1e-8;
        elseif local_tol==3
            thres=1e-10;
        end
        
        [fval,x,numeval]=dhc(fobj,X0,initsize,thres,100*nvar,x_L,x_U,weight,c_L,c_U,local_iterprint,tolc,varargin{:});
        exitflag=1;
        
    case 'ydhc'
        % adapted implementation of dhc algorithm
        
        nvar = numel(X0);
        
        if local_tol == 1
            thres = 1e-6;
        elseif local_tol == 2
            thres = 1e-8;
        else % local_tol == 3
            thres = 1e-10;
        end
        
        options.MaxFunEvals = 100*nvar;
        options.TolX = thres;
        options.TolFun = thres;
        
        if local_iterprint
            options.Display = 'iter';
        else
            options.Display = 'off';
        end
        
        fun = @(x) feval(fobj, x, varargin{:});
        
        [x, fval, exitflag, output] = ydhc(fun, X0, x_L, x_U, options);
        numeval = output.funcCount;
        
    case 'bobyqa'
        % matlab interface to Powell's bobyqa algorithm for
        % bound-constrained optimization
        
        nvar = numel(X0);
        
        options.MaxFunEvals = 100 * nvar;
        
        fun = @(x) feval(fobj, x, varargin{:});
        [x, fval, exitflag, output] = bobyqa(fun, X0, x_L, x_U, options);
        numeval = output.funcCount;
        
    case 'fsqp'
        %fsqp
        mode=110;
        
        iprint=0;   %Other values do not work
        
        nf=1;
        
        
        nineqn=length(n_upper)+length(n_lower);
        nineq=nineqn;
        
        neqn=neq;
        neq=neqn;
        
        
        miter=400;
        
        bigbnd=1e10;
        bl=x_L;
        bu=x_U;
        
        x0=X0';
        
        udelta=0;
        
        if local_tol==1
            tol1=tolc*100;
        elseif local_tol==2
            tol1=tolc;
        elseif local_tol==3
            tol1=tolc/100;
        end
        
        tol2=tol1;
        
        if nineq==0 && neq==0
            constr_file='';
        else
            constr_file='constr_fsqp';
        end
        
        [x,fval,g,lambda,info] = fsqp(mode,iprint,'fobj_fsqp',constr_file,'','',...
            nf,nineqn,nineq,neqn,neq,miter,udelta,bigbnd,...
            bl,bu,x0,tol1,tol2);
        
        
        numeval=n_fun_eval;
        x=x';
        
        if info==1 |info==2 |info==5 |info==6 | info==7 |info==9
            exitflag=0;
        else
            exitflag=1;
            
        end
        
        
    case 'ipopt'
        global params
        
        params=varargin;
        
        
        % The starting point.
        x0  = X0;  % The starting point.
        lb  = x_L;  % Lower bound on the variables.
        ub  = x_U;  % Upper bound on the variables.
        
        
        
        if isempty(c_U) & not(neq);
            nlc_fun='';
            nlc_jac='';
            lbc = [];   % Lower bounds on the constraint functions.
            ubc = [];   % Upper bounds on the constraint functions.
            
        else
            nlc_fun=@ipopt_c;
            nlc_jac=@eval_j_ipopt;
            lbc = c_L;   % Lower bounds on the constraint functions.
            ubc = c_U;   % Upper bounds on the constraint functions.
        end
        
        if local_iterprint
            dsp=5;
        else
            dsp=2;
        end
        
        [status x] = ipopt(x0,lb,ub,lbc,ubc,@ipopt_f,...
            @eval_g_ipopt,nlc_fun,...
            nlc_jac,'',[],'',[],'hessian_approximation','limited-memory',...
            'mu_strategy','adaptive',...                 % Update strategy for barrier parameter ['adaptive','monotone']
            'print_user_options','no',...                % Print all options set by the user
            'print_level',dsp,...                         % Output verbosity level [2 <= 5 <= 12 ]
            'nlp_scaling_method','gradient-based',...       % select the technique used for scaling the NLP ['none','user-scaling','gradient-based','equilibration-based']
            'print_options_documentation','no',...       % Switch to print all algorithmic options ['yes','no']
            'acceptable_dual_inf_tol',tolc,...% "Acceptance" threshold for the dual infeasibility [0 < 1e+10 < +inf]
            'acceptable_constr_viol_tol',tolc,... % "Acceptance" threshold for the constraint violation [0 < 0.01 < +inf]
            'max_iter',100,...                   % Maximal number of iterations allowed (Integer, >=0)
            'ma27_pivtol',tolc,...            % Pivot tolerance for the linear solvers [0 < 1e-08 < 1]
            'tol',tolc,...                    % Desired error tolerance (Double Precision, >0d0)[0 < 1e-08 < +inf]
            'acceptable_tol',tolc,...         % "Acceptable" convergence tolerance (relative) (Double Precision, >0d0) [0 < 1e-06 < +inf]
            'mu_init',1e-8,...                           % Initial value for the barrier parameter     [0 < 0.1 < +inf]
            'mu_min',tolc,...                 % Minimum value for barrier parameter [0 < 1e-11 < +inf]
            'bound_push',tolc,...             % Desired minimum absolute distance from the initial point to bound [0 < 0.01 < +inf]
            'bound_frac',tolc,...
            'derivative_test', 'none',...     % enable derivative checker ['none','first-order','second-order']);               % desired minimum relative distance from the initial point to bound [0 < 0.01 <= 0.5]
            'check_derivatives_for_naninf','no');    %Indicates whether it is desired to check for Nan/Inf in derivative matrices ['no','yes']
        
        
        
        fval=feval(fobj,x,varargin{:});
        exitflag=1;
        numeval=n_fun_eval+1;
        
    case 'misqp'
        
        ncont=length(X0)-int_var-bin_var;
        nint=int_var;
        nbin=bin_var;
        n=nint+ncont+nbin;
        m=length(n_upper)+length(n_lower)+neq;
        
        
        xl=x_L;
        xu=x_U;
        x0=X0;
        
        switch local_tol
            case 1
                acc=tolc*100;
            case 2
                acc=tolc;
            case 3
                acc=tolc/100;
        end
        
        [x,fval,nfunc,ifail,g]=run_misqp(ncont,nint,nbin,m,xl,xu,x0,acc,@fobj_misqp,neq,fobj,varargin{:});
        
        
        if ~isnan(fval) & ~isinf(fval) & all(g>=-tolc)
            exitflag=1;
        else
            exitflag=-1;
        end
        
        %        switch ifail
        %            case 0
        %                exitflag=1;
        %            case 4
        %                exitflag=-1;
        %            otherwise
        %                if all(g>=tolc)
        %                    exitflag=1;
        %                else
        %                    exitflag=-1;
        %                end
        %        end
        numeval=n_fun_eval;
        
        
    case 'fminsearch'
        switch local_tol
            case 1
                tolx=tolc*100;
                tolf=tolc*100;
            case 2
                tolx=tolc;
                tolf=tolc;
            case 3
                tolx=tolc/100;
                tolf=tolc/100;
            case 4
                tolx=tolc/10000;
                tolf=tolc/10000;
        end
        
        if local_iterprint
            dsp='iter';
        else
            dsp='off';
        end
        
        nvar=length(X0);
        
        options=optimset('LargeScale','off','Display',dsp,'Tolx',tolx,'TolFun',tolf,'MaxFunEvals',200*nvar,'MaxIter',1000*nvar);
        
        %[x,fval,exitflag,OUTPUT]=fminsearch(@fmobj,X0,options,fobj,[]);
        
        %[x,fval,exitflag,OUTPUT]=fminsearchbnd(@fmobj,X0,x_L,x_U,options,fobj,[],varargin{:});
        [x,fval,exitflag,OUTPUT]=fminsearchbnd(@fminsobj,X0,x_L,x_U,options,fobj,c_L,c_U,weight,tolc,varargin{:});
        numeval=n_fun_eval;
        
        
    case 'constrnew'
        switch local_tol
            case 1
                options(2)=tolc*100;
                options(3)=tolc*100;
                options(4)=tolc*100;
            case 2
                options(2)=tolc*1000;
                options(3)=tolc*1000;
                options(4)=tolc*1000;
            case 3
                options(2)=tolc/100;
                options(3)=tolc/100;
                options(4)=tolc/100;
            case 4
                options(2)=tolc/10000;
                options(3)=tolc/10000;
                options(4)=tolc/10000;
        end
        
        options(2)=1e-4;
        options(3)=1e-4;
        options(4)=tolc;
        
        
        if local_iterprint
            options(1)=1;
        else
            options(1)=-1;
        end
        
        options(13)=neq;
        
        nvar=length(X0);
        
        
        [x,OPTIONS,exitflag]=constrnew(@constr_obj,X0,options,x_L,x_U,[],fobj,neq,varargin{:});
        fval=OPTIONS(8);
        numeval=n_fun_eval;
        
        
        
    case 'lsqnonlin'
        switch local_tol
            case 1
                tolx=tolc;
                tolf=tolc;
            case 2
                tolx=tolc/100;
                tolf=tolc/100;
            case 3
                tolx=tolc/1000;
                tolf=tolc/1000;
            case 4
                tolx=tolc/100000;
                tolf=tolc/100000;
        end
        
        if local_iterprint
            dsp='iter';
        else
            dsp='off';
        end
        
        nvar=length(X0);
        
        % This set of options also works well
        %options=optimset('LargeScale','off','LevenbergMarquardt','off','LineSearchType','cubicpoly','Display',dsp,'Tolx',tolx,'TolFun',tolf);%,'MaxFunEvals',100*nvar);
        
        
        
        %If the solution goes beyond the bounds, use this set of options
        %         options=optimset('Display',dsp,'Tolx',tolx,'TolFun',tolf);%,'MaxFunEvals',100*nvar);
        
        % [x,resnorm,residual] = lsqnonlin(@lsqnonlin_fobj,X0,[],[],options,fobj,varargin{:});
        
        % This set of options (with largescale=on) to use the bounds
        if isempty(fjac)
            options=optimset('LargeScale','on','Display',dsp,'Tolx',tolx,'TolFun',tolf,'MaxIter',200);
        else
            % Jacobian is on:
            options=optimset('Jacobian','on','LargeScale','on','Display',dsp,'Tolx',tolx,'TolFun',tolf,'MaxIter',200);
        end
        fval=feval(fobj,X0,varargin{:});
        if isfinite(fval)
            [x,resnorm,residual] = lsqnonlin(@lsqnonlin_fobj,X0,x_L,x_U,options,fobj,fjac,varargin{:});
        else
            x = nan(size(X0));
            resnorm = inf;
            residual = [];
        end
        exitflag=1;
        fval=feval(fobj,x,varargin{:});
        numeval=n_fun_eval;
        
    case 'hooke'
        
        [x, histout] = hooke(X0', @hjobj,200*length(X0),tolc,local_iterprint,varargin{:});
        x=x';
        fval=histout(end,2);
        
        numeval=n_fun_eval;
        %        numeval=histout(end,1);
        
        exitflag=1;
        
    case 'fmincon'
        
        if isempty(c_U)
            const_fun=[];
        else
            const_fun=@fmcon;
        end
        
        switch local_tol
            case 1
                tolx=tolc*100;
                tolf=tolc*100;
                tolg=tolc*100;
            case 2
                tolx=tolc;
                tolf=tolc;
                tolg=tolc;
            case 3
                tolx=tolc/100;
                tolf=tolc/100;
                tolg=tolc/100;
            case 4
                tolx=tolc/10000;
                tolf=tolc/10000;
                tolg=tolc/10000;
        end
        % tolf=1e-9;
        % tolx=1e-9;
        if local_iterprint
            dsp='iter';
        else
            dsp='off';
        end
        %DW
        if local_opts.use_gradient_for_finish == 0
            
            options=optimset('LargeScale','off','Display',dsp,'Tolx',tolx,'TolFun',tolf,'Tolcon',tolg,'MaxSQPIter',100*length(X0),'MaxFunEvals',200*nvar,'MaxIter',200*nvar);
            
            [x,fval,exitflag,OUTPUT]=fmincon(@fmobj,X0,[],[],[],[],x_L,x_U,const_fun,options,fobj,neq,varargin{:});
            
        else
            
            %DW: provide gradient information
            if local_opts.check_gradient_for_finish
                options=optimset('LargeScale','off','Display',dsp,'Tolx',tolx,'TolFun',tolf,'Tolcon',tolg,'MaxSQPIter', ...
                    100*length(X0),'MaxFunEvals',200*nvar,'MaxIter',200*nvar, ...
                    'GradObj', 'on', 'DerivativeCheck', 'on'); %, 'SpecifyConstraintGradient', true);
            else
                options=optimset('LargeScale','off','Display',dsp,'Tolx',tolx,'TolFun',tolf,'Tolcon',tolg,'MaxSQPIter', ...
                    100*length(X0),'MaxFunEvals',200*nvar,'MaxIter',200*nvar, ...
                    'GradObj', 'on'); %, 'SpecifyConstraintGradient', true);
            end
            
            [x,fval,exitflag,OUTPUT]=fmincon(@fmobjgrad,X0,[],[],[],[],x_L,x_U,const_fun,options,fobj,neq,varargin{:});
            
        end
        
        numeval=n_fun_eval;
        
    otherwise
        error(['Local solver ' local_solver ' not recognized.']);
        
end
% make sure x is a row vector.
x = x(:).';

%\-----------------------------------------------------/
% Definition of constraints for constrnew
function [f,g] = constr_obj(x,fun,neq,varargin)
global n_fun_eval n_upper n_lower ccll ccuu

[f,ggg] = feval(fun,x,varargin{:});

g=ggg(1:neq);

ninequ=length(n_upper);
for i=1:ninequ
    g=[g ggg(n_upper(i))-ccuu(n_upper(i))];
end

nineql=length(n_lower);
for j=1:nineql
    g=[g ccll(n_lower(j))-ggg(n_lower(j))];
end

n_fun_eval=n_fun_eval+1;
return
%\-----------------------------------------------------/



%\-----------------------------------------------------/
% Definition of objective for lsqnonlin
function [fx, J] = lsqnonlin_fobj(x,fobj,fjac,varargin)
global n_fun_eval

if nargout > 1
    [ff,g,R] = feval(fobj,x,varargin{:});
    [GradObj,GradConstr, JacRes] = feval(fjac,x,varargin{:});
    J = JacRes;
else
    [ff,g,R] = feval(fobj,x,varargin{:});
end

fx=R;

n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/



%\-----------------------------------------------------/
% Definition of objective for fminsearchbnd
function fp = fminsobj(x,fun,c_L,c_U,weight,tolc,varargin)
global n_fun_eval

if isempty(who('n_fun_eval'))
    n_fun_eval=1;
end

if ~isempty(c_L) || ~isempty(c_U)
    [f,c] = feval(fun,x,varargin{:});
else
    [f] = feval(fun,x,varargin{:});
    c=[];
end

pen=ssm_penalty_function(x,c,c_L,c_U,tolc);

fp=f+weight*pen;

n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/



%\-----------------------------------------------------/
% Definition of objective for fmincon
function f = fmobj(x,fun,neq,varargin)
global n_fun_eval

f = feval(fun,x,varargin{:});
if isempty(who('n_fun_eval'))
    n_fun_eval=1;
end
n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/

%\-----------------------------------------------------/
% DW: Definition of objective for fmincon with gradient
function [f,g] = fmobjgrad(x,fun,neq,varargin)
global n_fun_eval

[f, ~, g] = feval(fun,x,varargin{:});
n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/

%\-----------------------------------------------------/
% Definition of constraints for fmincon
function [c,ceq] = fmcon(x,fun,neq,varargin)
global n_fun_eval n_upper n_lower ccll ccuu

c=[];
ceq=[];
[f,ggg] = feval(fun,x,varargin{:});

ceq=ggg(1:neq);

ninequ=length(n_upper);
for i=1:ninequ
    c=[c ggg(n_upper(i))-ccuu(n_upper(i))];
end

nineql=length(n_lower);
for j=1:nineql
    c=[c ccll(n_lower(j))-ggg(n_lower(j))];
end

n_fun_eval=n_fun_eval+1;
return
%\-----------------------------------------------------/



%\-----------------------------------------------------/
% Definition of objective and constraints for SOLNP
function f = snpfun(x,par,fun,nconst)
global n_fun_eval

if nconst
    [fobj,ci] = feval(fun,x,par{:});
    ci = ci(:);
else
    [fobj] = feval(fun,x,par{:});
    ci=[];
end
f = [fobj;
    ci];
n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/


%\-----------------------------------------------------/
% Definition of objective and constraints for misqp
function [fx,cx]= fobj_misqp(x,neq,fobj,varargin)

global n_fun_eval ccll ccuu n_upper n_lower

c=[];

if ~isempty(ccll) || ~isempty(ccuu)
    [fx,ggg] = feval(fobj,x,varargin{:});
else
    [fx] = feval(fobj,x,varargin{:});
    ggg=[];
end
n_fun_eval = n_fun_eval + 1;


ninequ=length(n_upper);
for i=1:ninequ
    c=[c ccuu(n_upper(i))-ggg(n_upper(i))];
end

nineql=length(n_lower);
for j=1:nineql
    c=[c ggg(n_lower(j))-ccll(n_lower(j))];
end

cx=[ggg(1:neq) c];

return
%\-----------------------------------------------------/

%
% %\-----------------------------------------------------/
% % Definition of objective for ipopt
% function f= ipopt_f(x,Params)
% global n_fun_eval
% fun=str2func(Params.fun);
% f=feval(fun,x);
% n_fun_eval = n_fun_eval + 1;
% return
% %\-----------------------------------------------------/
%
%
%
% %\-----------------------------------------------------/
% % Definition of constraints for ipopt
% function ceq = ipopt_c(x,dummy)
% global n_fun_eval
% fun=str2func(dummy.fun);
% neqc=dummy.neq;
% nvar=dummy.nvar;
% loc=dummy.loc;
%
% n_upper=loc.n_upper;
% n_lower=loc.n_lower;
% ccll=loc.ccll;
% ccuu=loc.ccuu;
%
% [f,ggg]=feval(fun,x);
% ceq=[];
% ceq=[ceq ggg(1:neqc)];
% clow=[];
% cupp=[];
%
%
% for i=1:length(n_upper)
%     cupp=[cupp; ggg(n_upper(i))-ccuu(n_upper(i)-neqc)+x(nvar+i)];
% end
%
% for j=1:length(n_lower)
%     clow=[clow; ggg(n_lower(j))-ccll(n_lower(j)-neqc)-x(nvar+length(n_upper)+j)];
% end
%
%
% ceq=[ceq cupp clow];
%
% n_fun_eval=n_fun_eval+1;
% return
% %\-----------------------------------------------------/
%
