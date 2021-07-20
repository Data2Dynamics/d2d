clear mex;
clear all;
close all;

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='ex3';                                    
problem.x_L=[0 0 0 0 0 0];      
problem.x_U=[1 1 1 1 16 16];  
problem.neq=4;
problem.c_L=-inf;
problem.c_U=4;

opts.ndiverse=25;

opts.local.solver='solnp';
opts.local.iterprint=1;
opts.local.tol=2;
%=========================================================================
k1=0.09755988;
k3=0.0391908;
k2=0.99*k1;
k4=0.9*k3;

Results_multistart=MEIGO(problem,opts,'MULTISTART',k1,k2,k3,k4);




