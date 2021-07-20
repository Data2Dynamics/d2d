% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssmgo_test.m 770 2013-08-06 09:41:45Z attila $
function ssmgo_test(solver,nproblem,noptim,pnames,lb,ub,prob,options,test,param);

% Function   : testssm
% Written by : Process Engineering Group IIM-CSIC (jegea@iim.csic.es)
% Created on : 20/04/2007
% Last Update: 21/10/2010
%
% Script to perform several optimizations with SSm/eSS for different
% problems 
% using the same or different settings
%
%          Calling:
%          ssmgo_test(nproblem,noptim,pnames,lb,ub,prob,options,test,param);
%
% Input Parameters:
%        
%                nproblem       = Number of problems to be tested (if you are testing the same problem n times, set nproblem=n, not 1).      
%                noptim         = Number of optimizations to perform per problem 
%                pames          = Cell array containing the names of the problems as strings
%                lb             = Cell array containing the lower bounds for all problems. 
%                ub             = Cell array containing the upper bounds for all problems. 
%                problem        = Matrix and structure containing problem settings for each problem.
%                                These settings are declared like in SSm but using indexes.
%                                For example, if we want to set an initial point for problem 3 we would type problem(3).x_0=[0 1 2 3 4];% 
%                opts           = General options for all problems. They are declared exactly the same way as in SSm
%                test           = Matrix and structure declaring specific options for problems.
%                param          = Structure declaring extra input parameters to be passed for every problem
% 
% Example: if the local solver chosen for all problems is fsqp except for problem 4 (for which we will choose misqp) we would write:
% 
%                       opts.local.solver='fsqp';
%                       test(4).local.solver='misqp';
%
% If problem 5 has 3 extra input arguments, p1, p2 and p3, we would declare
%       
%                       param{5}={p1,p2,p3};          
% 
% 
%  Output:
%
%   testssmgo generates two .mat files, one called Results_testssmgo_XXX.mat and testsummary_XXX.dat
%   (where XXX is a number related to the date and time when the test was performed), containing the following variables:
% 
%       Results_testssmgo_XXX.mat: Problem setting, options and Results (with the same out-puts as in SSm) 
%                                for each run under the format prob_px, opts_px and res_px_ry respectively,
%                                where x is the problem number and y is the run number.
% 
%       testsummary_XXX.mat: Summary of some results:
%                            
%                           fbest_px: Vector containing the best value found in each run for problem x.
%                           neval_px: Vector containing the number of function evaluations in each run for problem x.
%                           time_px: Vector containing the CPU time consumed in each run for problem x.
%                           best_values: Vector containing the best function value found for each problem after all the runs.
%                           worst_values: Vector containing the worst function value found for each problem after all the runs.
%                           mean_values: Vector containing the mean function value found for each problem after all the runs.
%
%
%

if nargin<10
    param=[];
    if nargin<9
        test=[];
    end
else
    aaa=length(param);
    if aaa<nproblem
        for i=aaa+1:nproblem
            param{i}=[];
        end
    end
end

n_times=0;
fff=[];
neval=[];
cpu_time=[];


hhh=clock;

flag=[strcat('Results_testssmgo_',num2str(hhh(1)),num2str(hhh(2)),num2str(hhh(3)),num2str(hhh(4)),num2str(hhh(5)))];
flag2=[strcat('testsummary_',num2str(hhh(1)),num2str(hhh(2)),num2str(hhh(3)),num2str(hhh(4)),num2str(hhh(5)))];

save(flag,'hhh');
save(flag2,'hhh');


fid=fopen('testssmgo_report.txt','w');
fprintf(fid,'%50s \n','SSMGO TEST FILE REPORT');
fprintf(fid,'%50s \n','**********************');
fprintf(fid,'%s \n',strcat(date));
fprintf(fid,'%s \n',  strcat(num2str(hhh(4)),':',num2str(hhh(5)),':',num2str(fix(hhh(6)))));
fprintf(fid,'Solver: %s \n',upper(solver));
fclose(fid);

best_values=[];
worst_values=[];
mean_values=[];

best_time=[];
worst_time=[];
mean_time=[];

best_neval=[];
worst_neval=[];
mean_neval=[];

size_test=length(test);

problems_solved=zeros(nproblem,noptim);

%Number of examples
for nfunc=1:nproblem
    clear mex
    tiempo=[];
    bestf=[];
    evals=[];

    
    problem=[];
    problem=prob(nfunc);
    problem.f=pnames{nfunc};
    problem.x_L=lb{nfunc};
    problem.x_U=ub{nfunc};


    opts=[];
    opts.local=[];
    opts=options;
    
    if not(isempty(test)) & nfunc<=size_test
        
        test_names=fieldnames(test(nfunc));

        if isfield(test(nfunc),'local')
            if not(isempty(test(nfunc).local))
                test_names_local=fieldnames(test(nfunc).local);
                for i=1:length(test_names_local)
                    opts.local.(test_names_local{i,:}) = test(nfunc).local.(test_names_local{i,:});
                end
            end
        end
        
        for i=1:length(test_names)
            if ~isempty(test(nfunc).(test_names{i,:}))
                opts.(test_names{i,:}) = test(nfunc).(test_names{i,:});
            end
        end
    end
    eval(['prob_p' num2str(nfunc) ' = problem;'])
    eval(['opts_p' num2str(nfunc) ' = opts;'])
    
    eval(['save(flag, [''prob_p'' num2str(nfunc)], ''-append'');'])
    eval(['save(flag, [''opts_p'' num2str(nfunc)], ''-append'');'])
    
    
    %Number of optimizations per example
    for j=1:noptim
        %Set seed for random numbers
        rand('state',j-1)

        if isempty(param) | isempty(param{nfunc})
            if strmatch(lower(solver),'ssm')
                Results=ssm_kernel(problem,opts);
            else
                Results=ess_kernel(problem,opts);
            end

        else
            if strmatch(lower(solver),'ssm')
                Results=ssm_kernel(problem,opts,param{nfunc}{:});
            else
                Results=ess_kernel(problem,opts,param{nfunc}{:});
            end
        end

        tiempo=[tiempo;Results.cpu_time];
        bestf=[bestf;Results.fbest];
        evals=[evals;Results.numeval];

        if Results.end_crit==3
            n_times=n_times+1;
            problems_solved(nfunc,j)=1;
        end

        fff=[fff Results.fbest];
        neval=[neval Results.numeval];
        cpu_time=[cpu_time Results.cpu_time];

        eval(['res_p' num2str(nfunc) '_r' num2str(j) ' = Results;'])
        eval(['save(flag, [''res_p'' num2str(nfunc) ''_r'' num2str(j)], ''-append'');'])
    end
    
    mean_tiempo=mean(tiempo);
    mean_bestf=mean(bestf);
    mean_evals=mean(evals);
    
    min_tiempo=min(tiempo);
    min_bestf=min(bestf);
    min_evals=min(evals);
    
    max_tiempo=max(tiempo);
    max_bestf=max(bestf);
    max_evals=max(evals);
    
    best_values=[best_values; min_bestf];
    worst_values=[worst_values; max_bestf];
    mean_values=[mean_values; mean_bestf];
    
    
    best_time=[best_time; min_tiempo];
    worst_time=[worst_time; max_tiempo];
    mean_time=[mean_time; mean_tiempo];
    
    
    best_neval=[best_neval; min_evals];
    worst_neval=[worst_neval; max_evals];
    mean_neval=[mean_neval; mean_evals];

    
    eval(['time_p' num2str(nfunc) ' = tiempo;'])
    eval(['fbest_p' num2str(nfunc) ' = bestf;'])
    eval(['neval_p' num2str(nfunc) ' = evals;'])
    
%     eval(['mean_time_p' num2str(nfunc) ' = mean_tiempo;'])
%     eval(['mean_fbest_p' num2str(nfunc) ' = mean_bestf;'])
%     eval(['mean_neval_p' num2str(nfunc) ' = mean_evals;'])
%     
%     eval(['min_time_p' num2str(nfunc) ' = min_tiempo;'])
%     eval(['min_fbest_p' num2str(nfunc) ' = min_bestf;'])
%     eval(['min_neval_p' num2str(nfunc) ' = min_evals;'])
%     
%     eval(['max_time_p' num2str(nfunc) ' = max_tiempo;'])
%     eval(['max_fbest_p' num2str(nfunc) ' = max_bestf;'])
%     eval(['max_neval_p' num2str(nfunc) ' = max_evals;'])


  
    eval(['save(flag2, [''time_p'' num2str(nfunc)],[''fbest_p'' num2str(nfunc)],[''neval_p'' num2str(nfunc)],''-append'');'])
    eval(['save(flag2, ''best_values'' ,''worst_values'',''mean_values'',''problems_solved'',''-append'');'])
    eval(['save(flag2, ''best_time'' ,''worst_time'',''mean_time'',''-append'');'])
    eval(['save(flag2, ''best_neval'' ,''worst_neval'',''mean_neval'',''-append'');'])
    eval(['save(flag2, ''solver'',''-append'');'])
%     eval(['save(flag2, [''mean_time_p'' num2str(nfunc)],[''mean_fbest_p'' num2str(nfunc)],[''mean_neval_p'' num2str(nfunc)],''-append'');'])
%     eval(['save(flag2, [''min_time_p'' num2str(nfunc)],[''min_fbest_p'' num2str(nfunc)],[''min_neval_p'' num2str(nfunc)],''-append'');'])
%     eval(['save(flag2, [''max_time_p'' num2str(nfunc)],[''max_fbest_p'' num2str(nfunc)],[''max_neval_p'' num2str(nfunc)],''-append'');'])
 

    

    fid=fopen('testssmgo_report.txt','a');
    fprintf(fid,'\n\n');
    fprintf(fid,'=========================================================================================================== \n');
    fprintf(fid,'%s %g \n','Problem',nfunc);
    fprintf(fid,'********** \n');
    fprintf(fid,'%s %s \n','Name: ',problem.f);
    fprintf(fid,'%s %g \n','nvar: ',length(problem.x_L));
    fprintf(fid,'\n');
    fprintf(fid,'RESULTS in %g optimizations \n',noptim);
    fprintf(fid,'*************************** \n');
    if isfield(problem,'vtr') & not(isempty(problem.vtr))
        fprintf(fid,'%s %g \n','Times the GO was found: ',n_times);
    end
    fprintf(fid,'%-10s  %-30s  %-50s  %-70s  \n','Run','Result','Nfuneval','CPU Time' );
    

    for k=1:length(fff)
        fprintf(fid,'%-10g  %-30e  %-50g  %-70g  \n',k,fff(k),neval(k)',cpu_time(k) );         
    end

    fprintf(fid,'\n\n');
    fprintf(fid,'%-25s %-40e  \n','Minimum function value:',min_bestf);
    fprintf(fid,'%-25s %-40e  \n','Maximum function value:',max_bestf);
    fprintf(fid,'%-25s %-40e  \n\n','Mean function value:',mean_bestf); 
    
    fprintf(fid,'%-25s %-40g  \n','Minimum func. eval.:',min_evals);
    fprintf(fid,'%-25s %-40g  \n','Maximum func. eval.:',max_evals);
    fprintf(fid,'%-25s %-40g  \n\n','Mean func. eval.:',round(mean_evals));
    
    fprintf(fid,'%-25s %-40g  \n','Minimum CPU time:',min_tiempo);
    fprintf(fid,'%-25s %-40g  \n','Maximum CPU time:',max_tiempo);
    fprintf(fid,'%-25s %-40g  \n\n','Mean CPU time:',mean_tiempo);
    
    fclose(fid);
    
    n_times=0;
    fff=[];
    neval=[];
    cpu_time=[];
    
end

fid=fopen('testssmgo_report.txt','a');

fprintf(fid,'************************************************************************************ \n');
fprintf(fid,'************************************************************************************ \n');

fprintf(fid,'%s \n','Number of solved problems');

for i=1:noptim
    aaa=find(problems_solved(:,i));
    fprintf(fid,'%s %i %s %i \n','Run',i,':',length(aaa));
end
fprintf(fid,'************************************************************************************ \n');
fprintf(fid,'************************************************************************************ \n');

bbb=sum(problems_solved');
n_diff_prob=numel(find(bbb));
fprintf(fid,'%s %i \n','Number of DIFFERENT solved problems: ',n_diff_prob);

total_solved_prob=sum(sum(problems_solved));
fprintf(fid,'%s %i \n','Number of TOTAL solved problems: ',total_solved_prob);

ccc=sum(problems_solved);
ddd=max(ccc);
eee=numel(find(ccc==ddd));

fprintf(fid,'%s %i %s %i %s \n','Maximum number of problems solved in a single run: ',ddd,' obtained ',eee,' times');


fclose(fid);
