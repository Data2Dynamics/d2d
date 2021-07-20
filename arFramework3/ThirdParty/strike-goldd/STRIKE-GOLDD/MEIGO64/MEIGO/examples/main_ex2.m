clear mex;
clear all;
close all;

% global optimum
%
%x*=[2.32952, 3.17849];
%    
%f(x*)=-5.50801


%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='ex2';                           %mfile containing the objective function
problem.x_L=[0 0];                         %lower bounds         
problem.x_U=[3 4];                         %upper bounds
problem.c_L=[-inf -inf];
problem.c_U=[2 36];

opts.maxeval=750;
opts.local.n1=2;
opts.local.n2=3;
%========================= END OF PROBLEM SPECIFICATIONS =====================
Results=MEIGO(problem,opts,'ESS');