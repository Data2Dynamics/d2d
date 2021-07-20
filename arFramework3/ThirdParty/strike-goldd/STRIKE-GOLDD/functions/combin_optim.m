%--------------------------------------------------------------------------
% Function that performs combinatorial optimization using the Variable
% Neighbourhood Search metaheuristic (VNS). It is a modification of the
% function 'rvnds_hamming' included in the MEIGO toolbox, introducing a
% time limit. Original 'rvnds_hamming' function created by Jose A. Egea:
% Egea, Jose A., et al. "MEIGO: an open-source software suite based on
% metaheuristics for global optimization in systems biology and 
% bioinformatics." BMC bioinformatics 15.1 (2014): 136.
%--------------------------------------------------------------------------

function [Results]=combin_optim(problem,opts)

%Max evaluations, time
maxeval=1e12;
maxtime=opts.maxtime;

%Initialize time
cpu_time=cputime;

%Load default values for the options
default=vns_defaults;

%Set all options
opts=vns_optset(default,opts);

%Set all options
opts=vns_optset(default,opts);
    
maxdist=opts.maxdist;
use_local=opts.use_local;
aggr=opts.aggr;
local_search_type=opts.local_search_type;
decomp=opts.decomp;
 
fobj = problem.f;
x_U  = problem.x_U;
x_L  = problem.x_L;

%Number of decision variables
nvar = length(x_L);

if not(isfield(problem,'x_0')) 
    x0=round(rand(1,nvar).*(x_U-x_L+1)+(x_L-0.5)); 
else
    x0=problem.x_0;
end

xcurrent=x0;
xbest=x0;

fbest=feval(fobj,x0);
nfuneval=1;
fcurrent=fbest;

fprintf('Init. point: Bestf: %g      CPUTime: %f\n',fbest,cputime-cpu_time);

Results.fbest=fbest;
Results.xbest=x0;
Results.func=fbest;
Results.time=cputime-cpu_time;
Results.neval=1;
Results.numeval=1;
Results.x=xbest;

%Initial neighborhood is k=1
kcurrent=1;

%Maximum neighborhood
kmax=floor(maxdist*nvar);

 %Perturb at least one dimension
 if not(kmax), kmax=1; end

improve=0;

while 1
     xnew=xcurrent;
     rand_var=randperm(nvar);
     shaked_vars=rand_var(1:kcurrent);
    for i=1:kcurrent
        continuar=0;
        tic
        while ~continuar && toc < 5 % limit time to avoid getting stuck
            xnew(shaked_vars(i))=round(rand*(x_U(shaked_vars(i))-...
                x_L(shaked_vars(i))+1)+(x_L(shaked_vars(i))-0.5));
            if xnew(shaked_vars(i))~=xcurrent(shaked_vars(i))
                continuar=1;
            end
        end
    end
    
    fnew=feval(fobj,xnew);
    
    nfuneval=nfuneval+1;
    
    
    if fnew<fbest
        fbest=fnew;
        xbest=xnew;
        xcurrent=xnew;
        fcurrent=fnew;
        improve=1;
        Results.func=[Results.func; fbest];
        Results.time=[Results.time; cputime-cpu_time];
        Results.neval=[Results.neval; nfuneval];
		Results.x=[Results.x;xbest];
    end
    
    if (cputime-cpu_time)>=maxtime
        
		fprintf('NFunEvals: %i  Bestf: %g      CPUTime: %f   \n',...
            nfuneval,fbest,cputime-cpu_time);
        fprintf('************************* \n')
        fprintf('END OF THE OPTIMIZATION \n')
        fprintf('Best solution value\t\t%g\n', fbest);
        fprintf('Decision vector\n');
        fprintf('\t%g\n', xbest');
        fprintf('CPU time\t\t%g\n', cputime-cpu_time);
        fprintf('Number of function evaluations\t\t%g\n\n', nfuneval);
        
        Results.xbest=xbest;
        Results.fbest=fbest;
        Results.func=[Results.func; fbest];
        Results.time=[Results.time; cputime-cpu_time];
        Results.neval=[Results.neval; nfuneval];
		Results.x=[Results.x;xbest];
        Results.numeval=nfuneval;
        Results.cpu_time=cputime-cpu_time;
        
		% save Results in a file  
		save VNS_report.mat Results problem opts  
        return
    end
       
    %Start the local phase 
    if use_local
        if ~aggr
		%In the non-aggressive scheme we apply local search over the new point even if it does no outperformed xbest
            f0=fnew;
            x0=xnew;
        elseif aggr && improve
            f0=fcurrent;
            x0=xcurrent;
        end
        
        if ~aggr || improve
		%In the aggressive scheme, it does not make sense to use the local search if xnew did not improve xbest
            % tic
             % fprintf('\nLOCAL SEARCH: Initial point function value: %g \n',...
                 % f0);
             [xout,fout,improve_local,evals_local,Results]=rvnds_local(x0,f0,...
                 fobj,x_L,x_U,local_search_type,shaked_vars,decomp,nvar,...
                 maxeval,maxtime,cpu_time,fbest,xbest,nfuneval, Results);
             % fprintf('Local solution function value: %g \n', fout);
             % fprintf('Number of function evaluations: %i \n', evals_local);
             % fprintf('Computation time: %f \n\n', toc);
             
            nfuneval=nfuneval+evals_local;
            if improve_local==2 %This means that the local search improved 
                                %the best value
                improve=1;
                xbest=xout;
                fbest=fout;
            end
        end
    end
    
    
    if improve
        improve=0;
        xcurrent=xbest;
        fcurrent=fbest;
        kcurrent=1;
        
        if ~use_local || improve_local<2
            fprintf('NFunEvals: %i  Bestf: %g      CPUTime: %f   \n',...
                nfuneval,fbest,cputime-cpu_time);
        end
        
    else
        kcurrent=kcurrent+1;
        if kcurrent>kmax
            kcurrent=1;
        end
    end
end
