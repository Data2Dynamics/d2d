% set parameter value by label
%
% arSetParsPattern(pattern, ... see arSetPars.m)
% 
% pattern	pattern of the parameter name

function arSetParsPattern(varargin)

global ar

pLabels = ar.pLabel(~cellfun(@isempty, strfind(ar.pLabel, varargin{1})));

arSetPars(pLabels, varargin{2:end});