% Model from J.W. Bolie. "Coefficients of normal blood glucose regulation". 
% J. Appl. Physiol., 16(5):783-788, 1961.

clear;

% 2 states
syms q1 q2
x = [q1; q2];

% 1 input
syms delta;
u = delta;

% 5 parameters 
syms p1 p2 p3 p4 Vp ic1 ic2
p =[p1; p2; p3; p4; Vp];

% 1 output
h = q1/Vp;

% initial conditions
ics  = [1,0]; 

% which initial conditions are known:
known_ics = [1,1];

% dynamic equations
f = [p1*q1-p2*q2+u;
    p3*q2+p4*q1];

save('Bolie_A_ics2','x','p','u','h','f','ics','known_ics');