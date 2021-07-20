clear mex;
clear all;
close all;


% global optimum
%
%x*=[0.0898, -0.7127];
% or    
%x*=[-0.0898, 0.7127];
%
%f(x*)= -1.03163;


	
	

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='ex1';                               %mfile containing the objective function
problem.x_L=-1*ones(1,2);                      %lower bounds         
problem.x_U=ones(1,2);                         %upper bounds

opts.maxeval=500;
opts.ndiverse=40;
opts.local.use_gradient_for_finish = 1; %DW: provide gradient to fmincon
opts.local.check_gradient_for_finish = 1; %DW: gradient checker
opts.local.solver='fmincon';
opts.local.finish='fmincon';
opts.local.iterprint=1;
%========================= END OF PROBLEM SPECIFICATIONS =====================

Results=MEIGO(problem,opts,'ESS');