function Initialize_FitTransient
global ar

ar.fit_transient.indp.signum    = strmatch('signum',ar.pLabel,'exact');
ar.fit_transient.indp.amp_sust  = strmatch('amp_sust',ar.pLabel,'exact');
ar.fit_transient.indp.amp_trans = strmatch('amp_trans',ar.pLabel,'exact');
ar.fit_transient.indp.offset    = strmatch('offset',ar.pLabel,'exact');
ar.fit_transient.indp.timescale_sust = strmatch('timescale_sust',ar.pLabel,'exact');
ar.fit_transient.indp.timescale_trans = strmatch('timescale_trans',ar.pLabel,'exact');
ar.fit_transient.indp.sd = strmatch('sd_Signal_au',ar.pLabel,'exact');

notSignum = setdiff(1:length(ar.pLabel),ar.fit_transient.indp.signum);
%  = setdiff(1:length(ar.pLabel),ar.fit_transient.indp.signum);
% 
% 

ar.ub(intersect(find(ar.qLog10),notSignum)) = 8;
ar.ub(intersect(find(~ar.qLog10),notSignum)) = 1000;


ar.qFit(ar.fit_transient.indp.signum) = 2;
ar.qLog10(ar.fit_transient.indp.offset) = 0;

ar.model.qPlotYs(:) = 0;
ar.model.qPlotYs(1) = 1;
