%--------------------------------------------------------------------------
% NF-kB transcription factor.
% The model is taken from:
%--------------------------------------------------------------------------
% Benjamin Merkt et al. (2015) Higher-order Lie symmetries in 
% identifiability and predictability analysis of dynamic models
%--------------------------------------------------------------------------

syms s1 s2 s3 s4...
    rhovol...
    k0 k1 k1p k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 x1_0 I0cyt I0nuc...
    x1 x2 x3 x4 x5 x6 x7 x8 x9 x10...
    y1 y2 y3 y4...
    u1...

% states:
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10].';

% outputs:
y1 = s1*(x1+x2+x3)+I0cyt;
y2 = s2*(x10+x5+x6)+I0nuc;
y3 = s3*(x2+x3);
y4 = s4*(x2+x4);
h = [y1 y2 y3 y4 ].';

% parameters:
p = [s1 s2 s3 s4 k0 k1 k1p k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 rhovol I0cyt I0nuc].';

% input:
w = u1;

% dynamic equations:
f = [k11*x10-(k1*u1/(1+k0*u1)+k1p)*x1;
    (k1*u1/(1+k0*u1)+k1p)*x1-k2*x2;
    k2*x2-k3*x3;
    k2*x2-k4*x4;
    k3*rhovol*x3-k5*x5;
    k5*x5-k10*x9*x6;
    k6*x6-k7*x7;
    k8*x7-k9*x8;
    k9*rhovol*x8-k10*x9*x6;
    k10*x9*x6-k11*rhovol*x10
    ];

% initial conditions:
ics  = [];
%ics = [x1_0,k1p*x1/k2,k1p*x1/k3,k1p*x1/k4,k3*rhovol*x3/k5,k7*x7/k6,k9*x8/k8,k3*x3/k9,k5*x5/(k10*x6),k1p*x1/k11];

% which initial conditions are known:
%known_ics = [1,0,0,0,0,0,0,0,0,1]; 
known_ics=[0,0,0,0,0,0,0,0,0,0];

save('NFKB_no_ics','x','p','h','f','w','known_ics');