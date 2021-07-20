clear all

%========================= PROBLEM SPECIFICATIONS ===========================
nvar=10;
problem.x_L=-5*ones(1,nvar);
problem.x_U=5*ones(1,nvar);
problem.f='rosen10';

opts.maxeval=1000;
algorithm='VNS';

[Results]=MEIGO(problem,opts,algorithm);
