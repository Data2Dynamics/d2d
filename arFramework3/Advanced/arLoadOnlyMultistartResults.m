% arLoadOnlyMultistartResults(workspace_name)
%
% Loads results of a multistart (e.g. LHS) fitting sequence
%
% Usage:
%       arLoadOnlyMultistartResults( workspace_name )
%
% Where:
%       workspace_name      Name of results folder
%
% See also arLoad

function arLoadOnlyMultistartResults(workspace_name)
if ~exist('Results','dir')
    error('No results folder exist. arLoad can only be executed in a D2D working directory.')
end
if exist('arCheck','file')==0
    arInit;
end
arCheck;

global ar


if(~exist('workspace_name', 'var') || isempty(workspace_name))
    [~, workspace_name] = fileChooser('./Results', 1, true);
elseif(isnumeric(workspace_name)) % workspace_name is the file-number
    [~, ~, file_list] = fileChooser('./Results', 1, -1);    
    workspace_name = file_list{workspace_name};
elseif(ischar(workspace_name)) 
    [~,workspace_name]=fileparts(workspace_name);    % remove path
end

tmpar = load(['./Results/' workspace_name '/workspace.mat']);

ar.psLabel = tmpar.ar.pLabel;%

ar.chi2s = tmpar.ar.chi2s;
ar.chi2sconstr = tmpar.ar.chi2sconstr;
ar.chi2s_start = tmpar.ar.chi2s_start;
ar.chi2sconstr_start = tmpar.ar.chi2sconstr_start;
ar.optim_crit = tmpar.ar.optim_crit;
ar.ps = tmpar.ar.ps;
ar.ps_start = tmpar.ar.ps_start;
ar.exitflag = tmpar.ar.exitflag;


% from arPlotFits
chi2constr = ar.chi2s + ar.chi2sconstr;
[~, isorted] = sort(chi2constr);
ar.chi2s_sorted = ar.chi2s(isorted);
ar.chi2s_start_sorted = ar.chi2s_start(isorted);
ar.chi2sconstr_sorted = ar.chi2sconstr(isorted);
ar.chi2sconstr_start_sorted = ar.chi2sconstr_start(isorted);
ar.ps_sorted = ar.ps(isorted,:);
ar.ps_start_sorted = ar.ps_start(isorted,:);


fprintf('LHS loaded and sorted from file %s\n', workspace_name);
end




