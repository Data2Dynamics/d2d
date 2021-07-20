% Model from J.W. Bolie. "Coefficients of normal blood glucose regulation". 
% J. Appl. Physiol., 16(5):783-788, 1961.
% Input function unknown
clear;
% 2 states
syms q1 q2
x = [q1; q2];

% 1 input
syms delta;
w = delta;

% 5 parameters 
syms p1 p2 p3 p4 Vp 
p =[p1; p2];

% 1 outputf
h = q1/Vp;

% initial conditions
ics  = []; 

% which initial conditions are known:
known_ics = [0,0];

% dynamic equations
f = [p1*q1-p2*q2+w;
    p3*q2+p4*q1];

save('Bolie_A_cha','x','p','w','h','f','ics','known_ics');