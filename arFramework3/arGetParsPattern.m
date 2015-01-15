% get parameter value by pattern
%
% [p, pLabels] = arGetParsPattern(pattern, qLog10)
% 
% pattern	pattern of the parameter name
% qLog10	0=normal, 1=log10 parameter values

function [p, pLabels] = arGetParsPattern(pattern, qLog10)

global ar

pLabels = ar.pLabel(~cellfun(@isempty, strfind(ar.pLabel, pattern)));

p = arGetPars(pLabels, qLog10);
