clear mex;
clear all;
close all;

%global optimum
%
%x*=[0.77152
% 	0.516994
% 	0.204189
% 	0.388811
% 	3.0355
% 	5.0973];
% 
%f(x*)= -0.388811;


%====================== PROBLEM DEPENDENT DECLARATIONS=====================
%=================== END OF PROBLEM DEPENDENT DECLARATIONS ================


%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='ex3';                                    
problem.x_L=[0 0 0 0 0 0];      
problem.x_U=[1 1 1 1 16 16];  
problem.neq=4;
problem.c_L=-inf;
problem.c_U=4;
   
opts.maxtime=5;
opts.local.solver='solnp';
%========================= END OF PROBLEM SPECIFICATIONS =====================
k1=0.09755988;
k3=0.0391908;
k2=0.99*k1;
k4=0.9*k3;

[Results]=MEIGO(problem,opts,'ESS',k1,k2,k3,k4);
