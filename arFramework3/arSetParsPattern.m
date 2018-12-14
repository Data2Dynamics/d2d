% arSetParsPattern(pattern, options)
% 
% pattern	pattern of the parameter name
% options   works in the same manner as arSetPars. See there for more
%           information
% 
% see also arSetPars

function arSetParsPattern(varargin)

global ar

pLabels = ar.pLabel(~cellfun(@isempty, strfind(ar.pLabel, varargin{1})));

arSetPars(pLabels, varargin{2:end});