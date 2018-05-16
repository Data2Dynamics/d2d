function[model]=AktPathway_syms()

% AktPathway_syms for examples/aktPathway
%
% creates an amimodel-object for the AMICI solver
%
% Parameters: 
%  void
% 
% Return values:
%  model: amimodel object


%% STATES
syms pro_EGFR EGFR EGF_EGFR pEGFR Akt pEGFR_Akt pAkt S6 pAkt_S6 pS6 EGF
x=[pro_EGFR EGFR EGF_EGFR pEGFR Akt pEGFR_Akt pAkt S6 pAkt_S6 pS6 EGF];

%% PARAMETERS:
syms re3_k1 re3_k2 re4_k1 re5_k1 re6_1_k1 re6_1_k2 re6_2_k1 re7_k1 re8_1_k1 re8_1_k2 re8_2_k1 re9_k1 EGFR_turnover pEGFR_scaleFactor pAkt_scaleFactor pS6_scaleFactor;
p=[re3_k1 re3_k2 re4_k1 re5_k1 re6_1_k1 re6_1_k2 re6_2_k1 re7_k1 re8_1_k1 re8_1_k2 re8_2_k1 re9_k1 EGFR_turnover pEGFR_scaleFactor pAkt_scaleFactor pS6_scaleFactor];

%% CONDITIONS:
syms EGF_conc_step EGF_conc_impulse EGF_conc_ramp pulse_time ramp_time;
k=[EGF_conc_step EGF_conc_impulse EGF_conc_ramp pulse_time ramp_time];

%% DYNAMICS:
syms t
f(2)=re3_k2*EGF_EGFR - re3_k1*EGF*EGFR + EGFR_turnover *(pro_EGFR - EGFR);
f(3)=re3_k1*EGF*EGFR - re3_k2*EGF_EGFR - re4_k1*EGF_EGFR;
f(4)=re4_k1*EGF_EGFR - re5_k1*pEGFR + re6_2_k1*pEGFR_Akt - re6_1_k1*pEGFR*Akt + re6_1_k2*pEGFR_Akt;
f(6)=re6_1_k1*pEGFR*Akt - re6_1_k2*pEGFR_Akt - re6_2_k1*pEGFR_Akt;
f(5)=re7_k1*pAkt - re6_1_k1*pEGFR*Akt + re6_1_k2*pEGFR_Akt;
f(7)=re6_2_k1*pEGFR_Akt - re8_1_k1*pAkt*S6 + re8_1_k2*pAkt_S6 + re8_2_k1*pAkt_S6 - re7_k1*pAkt;
f(9)=re8_1_k1*pAkt*S6 - re8_1_k2*pAkt_S6 - re8_2_k1*pAkt_S6;
f(8)=re9_k1*pS6 - re8_1_k1*pAkt*S6 + re8_1_k2*pAkt_S6;
f(10)=re8_2_k1*pAkt_S6 - re9_k1*pS6 ;
f(11)=EGF_conc_ramp/ramp_time; %gives warning for ramp_time=0, but replaces value with 0.0
f(1)=0;


%% INITIAL CONDITIONS:
x0=[68190.183733379701,68190.183733379701,0,0,0.043309016570930899,0,0,3.5431674054221798,0,0,EGF_conc_step+EGF_conc_impulse];

%% EVENTS:
event=amievent(t-pulse_time,[0,0,0,0,0,0,0,0,0,0,-EGF_conc_impulse],[]);

%% OBSERVABLES:
%pEGFR_total:
y(1)=(pEGFR+pEGFR_Akt)*pEGFR_scaleFactor;
%pAkt_total
y(2)=(pAkt + pAkt_S6)*pAkt_scaleFactor;
%pS6_total:
y(3)=pS6*pS6_scaleFactor;



model.sym.p = p;
model.sym.k = k;
model.sym.x = x;
model.sym.y = y;
model.sym.xdot = f;
model.param = 'log10';
model.event = event;
model.sym.x0 = x0;

end
