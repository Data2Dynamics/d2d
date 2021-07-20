% Nonlinear model with dynamical compensation 
% Originally published in: Karin et al, Mol Syst Biol 2016
% Corresponds to the model in Fig. 1C of: Villaverde & Banga, arXiv:1701.02562

clear;

% 3 states
syms x1 x2 x3 
x = [x1; x2; x3];

% 1 output
h = x1;

% 1 input
syms u u0;

% 2 parameters 
syms p1 p2 
p =[p1; p2];

% initial conditions
syms x10 x20 x30
ics  = [x10 x20 x30];

% which initial conditions are known:
known_ics = [1,1,1];

% dynamic equations
f = [u-p2*x1*x3;
    x2*(x1-x10);
    p1*x2*x1-x3];


save('1C_nonlinear_known_ics_par','x','p','u','h','f','ics','known_ics');