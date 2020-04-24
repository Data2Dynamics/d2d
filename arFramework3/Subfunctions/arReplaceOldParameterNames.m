% ar = arReplaceOldParameterNames(ar)
% 
%   This function can be edited in order to replace parameter names
%   (i.e. ar.pLabel). This function is called by D2D in arLoadPars and
%   pleCompare
% 
%   The following procedure has to be applied:
%   1) A local copy of this function has to be used (that has higher priority
%   in Matlab's search path). 
%   2) Replacements have to be done as specified below
% 
% See also arLoadPars, pleCompare

function pLabel = arReplaceOldParameterNames(pLabel)

% pLabel = strrep(pLabel,'OldName','newName');

