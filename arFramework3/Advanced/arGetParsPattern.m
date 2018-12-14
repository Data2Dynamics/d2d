% [p, pLabels] = arGetParsPattern(pattern, qLog10)
%
% get parameter values by matching pattern to ar.pLabel
% 
% pattern	parameter names are searched for this pattern
% qLog10	logical to get log10 of parameter value [false]
%
% p         parameter values with found pattern in it
% pLabels   cell array of parameter names with found pattern in it
%
% See also arGetPars arPrint

function [p, pLabels] = arGetParsPattern(pattern, qLog10)

global ar

pLabels = ar.pLabel(~cellfun(@isempty, strfind(ar.pLabel, pattern)));

p = arGetPars(pLabels, qLog10);
