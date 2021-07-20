%--------------------------------------------------------------------------
% The model is taken from:
% Raia V et al. (2011) Dynamic mathematical modeling of IL13-induced 
% signaling in Hodgkin and primary mediastinal B-cell lymphoma allows 
% prediction of therapeutic targets. 
% Cancer Research 71(3):693{704.
%--------------------------------------------------------------------------
% The model contains certain modifications in relation to the 
% initial equations
%--------------------------------------------------------------------------
clear all;

% 10 states:
syms x1 x2 x3 x4 x5 x6 x8 x10 x11 x13
x = [x1; x2; x3; x4; x5; x6; x8; x10; x11; x13];

% 23 parameters (t23 is an initial condition):
syms t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12...
    t13 t14 t15 t16 t17 t18 t19 t20 t21 t22 t23
p = [t1; t2 ; t3 ; t4 ; t5 ; t6 ; t7 ; t8 ; t9 ; t10 ; t11 ; t12...
    ; t13 ; t14 ; t15 ; t16 ; t17; t18 ; t19 ; t20 ; t21 ; t22];

% 8 outputs:
h = [
x1 + x3 + x4;    
t18*(x3 + x4 + x5 +(0.34-x11));
t19*(x4 + x5);
t20*(-x6+2.8);
t21*x10;
t22*x10*t17/t11;
x13;
-x8+165
];

% 2 known constants:
c1 = 2.265;
c2 = 91;

% one input:
syms u1;
u = u1;

% dynamic equations:
f = [ 
-t1*x1*c1*u1-t5*x1+t6*x2;
t5*x1-t6*x2;
t1*c1*u1*x1-t2*x3*(-x6+2.8);
t2*x3*(-x6+2.8)-t3*x4;
t3*x4-t4*x5;
-t7*x3*x6*t11/(t11+x10*t13*t17)-t7*x4*x6*t11/(t11+x10*t13*t17)+t8*(-x6+2.8)*c2;
-t9*x8*(-x6+2.8)+t10*(-x8+165)*c2;
t11*(-x8+165);
-t12*c1*u1*x11;
x10*t14/(t15+x10)-t16*x13
];

% initial conditions:
ics  = [1.3, t23, 0, 0, 0, 2.8, 165, 0, 0.34, 0];   
%ics=[];

% which initial conditions are known:
known_ics = [1,1,1,1,1,1,1,1,1,1]; 
%known_ics = [0,0,0,0,0,0,0,0,0,0]; 

save('JAKSTAT_modified','x','p','h','f','u','ics','known_ics');
