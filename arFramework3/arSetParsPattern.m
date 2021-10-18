% arSetParsPattern(pattern, options)
% 
% pattern	pattern of the parameter name
% {pattern, n} maximum number of characters to compare
% options   works in the same manner as arSetPars. See there for more
%           information
%
% Example:
% arSetParsPattern('a_',5,1,1,-5,-2)
% arSetParsPattern({'a_',2},5,1,1,-5,-2) 
% 2nd option: just if pattern is in the beginning of the parameter name
% 
% see also arSetPars

function arSetParsPattern(varargin)

global ar

if iscell(varargin{1})
    pLabels = ar.pLabel(strncmp(ar.pLabel, varargin{1}{1},varargin{1}{2}));
else
    pLabels = ar.pLabel(~cellfun(@isempty, strfind(ar.pLabel, varargin{1})));
end

arSetPars(pLabels, varargin{2:end});