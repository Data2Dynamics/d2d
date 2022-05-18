% Input functions, usually depending on the independent variable, can be specified.
% Inputs do not belong to a specific compartment. The user has to ensure plausibility 
% of the usage of the input function himself.
%
% Example: 
%   INPUTS
%   Epo             C   "units/cell"   "conc."   "k1*exp(-k2*t)"
%   SAv             C   "units/cell"   "conc."   "k3"
%
% The 1st argument specifies the unique identifier (here Epo and SAv) that can be 
% used in the mathematical expressions later.
% The 2nd-4th argument specify the unit types (here C = concentration), the unit 
% itself (here "units/cell"), and a plain text used for plots (here "conc.").
% The 5th argument specify the mathematical expression of the input function 
% (here "k1*exp(-k2*t)" is an exponential decay and "k3" is a constant with 
% unknown value).
%
% STEP FUNCTIONS
%   "step1(t, level1, switch_time, level2)" 
%     Provides a simple step function at the value switch_time of the independent variable (usually time).
%   "step2(t, level1, switch_time1, level2, switch_time2, level3)" 
%     Provides a double step function at the values switch_time1/2 of the independent variable t.
%   "step3(t, level1, switch_time1, level2, switch_time2, level3, switch_time3, level4)" 
%     Provides a double step function at the values switch_time1/2/3 of the independent variable t.
%   "step4(t, level1, switch_time1, level2, switch_time2, level3, switch_time3, level4, switch_time4, level5)" 
%     Provides a double step function at the values switch_time1/2/3/4 of the independent variable t.
%   "smoothstep1(t, level1, switch_time, level2, smoothness)" 
%     Provides a smooth step function at the value switch_time of the independent variable (usually time).
%   "smoothstep2(t, level1, switch_time1, level2, switch_time2, level3, smoothness)" 
%     Provides a smooth double step function at the values switch_time1/2 of the independent variable t.
%     The smoothness parameter controls the smoothness of the step. For large values, the step becomes increasingly sigmoidal.
%
%   *** Must use arFindInputs after compilation ***
% 
% BOLUS INJECTIONS
%   "bolus(t, amount, time_point, duration)" 
%     giving a bolus injection at value time point of the independent variable t. This corresponds to a Gaussian curve with area amount and standard deviation duration.
%
% SPLINES
%   Cubic splines with 3, 4, 5 or 10 knots can be defined by:
%     "spline3(t, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3, q_initial_slope_constraint, initial_slope)"
%     "spline4(t, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3, t_knot4, p_knot4, q_initial_slope_constraint, initial_slope)"
%     "spline5(t, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3, t_knot4, p_knot4, t_knot5, p_knot5, q_initial_slope_constraint, initial_slope)"
%     "spline10(t, t_knot1, p_knot1, t_knot2, p_knot2, ..., t_knot10, p_knot10, q_initial_slope_constraint, initial_slope)"
%       Here, t denotes the independent variable, t_knot1-n indicate the locations of the spline knots (these should be fixed to a numeric value), 
%       p_knot1-5 denote the spline parameters (these can either be fixed or left as free parameters to be estimated) and q_initial_slope_constraint 
%       is a flag (0 = no, 1 = yes) that indicates if the slope of the spline at the first knot is to be constrained to the value given in initial_slope 
%       (both values should be fixed to a numeric value). The spline parameter p_knot is defined as the value u(t_knot) of the spline at t_knot. If the 
%       spline should be constrained to positive values use the functions spline_pos3, spline_pos4 and spline_pos5 the same way as described above. 
%       An example using splines is given in Examples/Swameye_PNAS2003.
%
%   Monotonic splines with 3, 4, 5 or 10 knots can be defined by:
%     "monospline3(t, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3)"
%     "monospline4(t, t_knot1, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3, t_knot4, p_knot4)"
%     "monospline5(t, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3, t_knot4, p_knot4, t_knot5, p_knot5)"
%     "monospline10(t, t_knot1, p_knot1, t_knot2, p_knot2, t_knot3, p_knot3, ..., t_knot9, p_knot9, t_knot10, p_knot10)"
%       Monotonic are splines which are monotonic between knots. Their syntax is similar to that of the regular splines. 
%       Monotonic splines have as an extra advantage that they do not suffer from over and undershoot behaviour as much 
%       as the regular cubic splines, considering that the interpolating function between two knots is guaranteed to be 
%       monotonic. They are useful for modelling input sources which are not noisy, or inputs requiring a large number 
%       of knots. Setting the first and last two knot parameters to the same value guarantees that the monotonic spline 
%       remains constant outside the spline range. This means that monotonic splines can be combined into splines comprising 
%       of more than 10 knots (see Examples/LongSplines).
%
%   Input splines
%     "inputspline(t, N, [t_knot1, t_knot2, ... t_knotN], [p_knot1, p_knot2, ... p_knotN])"
%       The inputspline is a special version of the monotonic spline. It allows for a spline of arbitrary length, 
%       but *cannot be fitted* (no parameters are allowed in the time or parameter series).
%
%  Note: When the splines are slow, invoke ar.config.turboSplines = 1, to use a more performant spline implementation.
%
% GENERAL INPUT FUNCTIONS
%   Often, an input consists of transient and sustained parts. Such behaviour can be implemented by the following expression:
%     "gif_amp_trans*(1-exp(-t/gif_timescale_sust))*exp(-t/(gif_timescale_trans)) + gif_amp_sust*(1-exp(-t/gif_timescale_sust))"
%       The function has three parameters, two amplitudes (gif_amp_trans and gif_amp_sust) and two time scales (gif_timescale_trans 
%       and gif_timescale_sust), that encode the transient and sustained parts.

help arhelpInputs