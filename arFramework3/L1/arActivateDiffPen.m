function y = arActivateDiffPen(pairs)

% function y = arActivateDiffPen(pairs)
%
% Activates the symmetric penalization of fold-change differences. Use
% before calling arRegularize. The argument 'pairs' is an output of 
% arPrintDiffPenConditions.
%
%   pairs   Cell-array of parameter names with each row containing two
%           parameters of which the difference is penalized, as well as the
%           corresponding boolean isEq-Parameter, that toggles constraint
%           fitting on or off.
%
% 
% Example:  arActivateDiffPen({'relto_CL2_p1', 'relto_CL3_p1', 'isEq_CL2_CL3_p1';
%           'relto_CL2_p2', 'relto_CL3_p2', 'isEq_CL2_CL3_p2'});
%           arRegularize('relto',0,4,'Type','lq','Deform',1-0.8)
%
% See also: arPrintDiffPenConditions, arRegularize

global ar

N = size(pairs,1);
L1DiffPen_diffs = NaN(N,2);
L1DiffPen_constrSwitches = NaN(N,1);
for i = 1:N
    L1DiffPen_diffs(i,1) = arFindPar(pairs{i,1},'exact');
    L1DiffPen_diffs(i,2) = arFindPar(pairs{i,2},'exact');
    L1DiffPen_constrSwitches(i) = arFindPar(pairs{i,3}, 'exact');
end

ar.L1DiffPen_diffs = L1DiffPen_diffs;
ar.L1DiffPen_constrSwitches = L1DiffPen_constrSwitches;
ar.L1DiffPen_activate = true;
ar.L1DiffPen_useCustomRes = false;
end