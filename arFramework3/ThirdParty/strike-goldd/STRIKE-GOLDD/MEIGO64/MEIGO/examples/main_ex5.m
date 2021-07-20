clear mex;
clear all;
close all;


% global optimum
%
%x*=[5.92585e-005
%   2.9634e-005
%   2.0471e-005
%   0.000274465
%   3.99821e-005];
%    
%f(x*)=19.8721

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='ex5';           

problem.x_L=zeros(1,5);
problem.x_U=ones(1,5);

problem.x_0=0.5*ones(1,5);
opts.maxeval=1e3; 
opts.log_var=[1:5];
opts.local.solver='nl2sol';
opts.inter_save=1;
%========================= END OF PROBLEM SPECIFICATIONS =====================
%time intervals

texp=[0 1230 3060 4920 7800 10680 15030 22620 36420];
 
% Distribution of species concentration
 %	     y(1)    y(2)    y(3)    y(4)    y(5)
 
 yexp=[ 100.0	 0.0	 0.0	 0.0     0.0
    	88.35    7.3     2.3     0.4     1.75
        76.4    15.6     4.5     0.7     2.8
        65.1    23.1     5.3     1.1     5.8
        50.4    32.9     6.0     1.5     9.3
        37.5    42.7     6.0     1.9    12.0
        25.9    49.1     5.9     2.2    17.0
        14.0    57.4     5.1     2.6    21.0
         4.5    63.1     3.8     2.9    25.7 ];


Results=MEIGO(problem,opts,'ESS',texp,yexp);