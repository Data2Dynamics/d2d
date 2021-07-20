function [R]= objf_nl2sol(x)
% Objective function for NL2SOL in AMIGO
% Automatically generated in ssm_aux_local.m
global n_fun_eval 
global input_par 
[f,g,R]= ex5(x,input_par{:});
n_fun_eval=n_fun_eval+1; 
return
