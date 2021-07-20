function [Results]=rvnds_hamming(problem,opts,varargin);

if not(isfield(opts,'maxeval')) & not(isfield(opts,'maxtime'))
    fprintf('WARNING:Either opts.maxeval or opts.maxtime must be defined as a stop criterion \n')
    fprintf('Define any of these options and rerun \n')
    Results=[];
    return
else
    if not(isfield(opts,'maxeval'))
        maxeval=1e12;
    else
        maxeval=opts.maxeval;
    end
    if not(isfield(opts,'maxtime'))
        maxtime=1e12;
    else
        maxtime=opts.maxtime;
    end
end


%Initialize time
cpu_time=cputime;

%Load default values for the options
default=vns_defaults;

if nargin<2, opts=[]; end

%Set all options
opts=vns_optset(default,opts);


%Set all optiions
opts=vns_optset(default,opts);
    
maxdist=opts.maxdist;
use_local=opts.use_local;
aggr=opts.aggr;
local_search_type=opts.local_search_type;
decomp=opts.decomp;
 
fobj=problem.f;
x_U=problem.x_U;
x_L=problem.x_L;


%Check if bounds have the same dimension
if length(x_U)~=length(x_L)
    disp('Upper and lower bounds have different dimension!!!!')
    disp('EXITING')
    Results=[];
    return
else
    %Number of decision variables
    nvar=length(x_L);
end

if not(isfield(problem,'x_0')), x0=round(rand(1,nvar).*(x_U-x_L+1)+(x_L-0.5)); else, x0=problem.x_0; end

xcurrent=x0;
xbest=x0;

fbest=feval(fobj,x0,varargin{:});
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
        while ~continuar
            xnew(shaked_vars(i))=round(rand*(x_U(shaked_vars(i))-...
                x_L(shaked_vars(i))+1)+(x_L(shaked_vars(i))-0.5));
            if xnew(shaked_vars(i))~=xcurrent(shaked_vars(i));
                continuar=1;
            end
        end
    end
    
    fnew=feval(fobj,xnew,varargin{:});
    
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
    
    if nfuneval>=maxeval | (cputime-cpu_time)>=maxtime
        
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
        elseif aggr & improve
            f0=fcurrent;
            x0=xcurrent;
        end
        
        if ~aggr | improve
		%In the aggressive scheme, it does not make sense to use the local search if xnew did not improve xbest
            % tic
             % fprintf('\nLOCAL SEARCH: Initial point function value: %g \n',...
                 % f0);
             [xout,fout,improve_local,evals_local,Results]=rvnds_local(x0,f0,...
                 fobj,x_L,x_U,local_search_type,shaked_vars,decomp,nvar,...
                 maxeval,maxtime,cpu_time,fbest,xbest,nfuneval, Results, varargin{:});
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
        
        if ~use_local | improve_local<2
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
