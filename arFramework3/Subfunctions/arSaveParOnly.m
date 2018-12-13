%   arSaveParOnly(ar)
%
%   arSaveParOnly(ar, savepath, filename)
%
%
% This function is called by arSave or can be used independently of arSave.
% It extracts the important fields from the ar struct and saves it as
% workspace in savepath/workspace
% save only parameters
%
%   savepath        Default: ar.config.savepath
%
%   filename        Default: 'workspace_pars_only.mat'
%
%   Please note, that the global ar is not used here. Instead, ar is
%   provided as a variable because the variable ar is is required for
%   storing the compressed version.
%

function arSaveParOnly(arIn, savepath, filename)
if ~exist('savepath','var') || isempty(savepath)
    savepath = arIn.config.savepath;
end

if ~exist('filename','var') || isempty(filename)
    filename = 'workspace_pars_only.mat';
end

% check if dir exists
if(~exist(savepath, 'dir'))
    mkdir(savepath)
end

arIn = arUpdateCheckstr(arIn);

ar = struct([]);
ar(1).pLabel = arIn.pLabel;
ar.p = arIn.p;
ar.qLog10 = arIn.qLog10;
ar.qFit = arIn.qFit;
ar.lb = arIn.lb;
ar.ub = arIn.ub;
if isfield(arIn,'type') % not available in older versions
    ar.type = arIn.type;
end
if isfield(arIn,'mean') % not available in older versions
    ar.mean = arIn.mean;
end
if isfield(arIn,'std') % not available in older versions
    ar.std = arIn.std;
end
if isfield(arIn,'fkt') % not available in older versions
    ar.fkt = arIn.fkt;
end
if isfield(arIn,'checkstrs') % not available in older versions
    ar.checkstrs = arIn.checkstrs;
end

try %#ok<TRYNC>
    ar.chi2fit = arIn.chi2fit;
    ar.ndata = arIn.ndata;
    ar.nprior = arIn.nprior;
end
try %#ok<TRYNC>
    ar.ps = arIn.ps;
    ar.ps_errors = arIn.ps_errors;
    ar.chi2s = arIn.chi2s;
    ar.chi2sconstr = arIn.chi2sconstr;
    ar.timing = arIn.timing;
    ar.exitflag = arIn.exitflag;
    ar.fun_evals = arIn.fun_evals;
    ar.ps_start = arIn.ps_start;
    ar.chi2s_start = arIn.chi2s_start;
    ar.chi2sconstr_start = arIn.chi2sconstr_start;
    ar.optim_crit = arIn.optim_crit;
    
    ar.chi2s_sorted = arIn.chi2s_sorted;
    ar.chi2sconstr_sorted = arIn.chi2sconstr_sorted;
    ar.ps_sorted = arIn.ps_sorted;
    ar.chi2s_start_sorted = arIn.chi2s_start_sorted;
    ar.chi2sconstr_start_sorted = arIn.chi2sconstr_start_sorted;
    ar.ps_start_sorted = arIn.ps_start_sorted;
end
if isfield(arIn,'config')
    if isfield(arIn.config,'useFitErrorMatrix') && arIn.config.useFitErrorMatrix == 0
        ar.config.fiterrors = arIn.config.fiterrors; % #ok<STRNU>
        ar.config.useFitErrorCorrection = arIn.config.useFitErrorCorrection;
    elseif isfield(arIn.config,'fiterrors_matrix')
        ar.config.useFitErrorMatrix = true;
        ar.config.fiterrors_matrix = arIn.config.fiterrors_matrix;
        ar.config.ploterrors_matrix = arIn.config.ploterrors_matrix;
        ar.config.useFitErrorCorrection = arIn.config.useFitErrorCorrection;
    elseif isfield(arIn.config,'fiterrors')  % the usual case
        ar.config.fiterrors = arIn.config.fiterrors;
        if isfield(arIn.config,'useFitErrorCorrection')
            ar.config.useFitErrorCorrection = arIn.config.useFitErrorCorrection;
        end
    end
end

if(isfield(arIn,'ple'))
    if (isfield(arIn.ple,'chi2s'))
        ar.ple.chi2s = arIn.ple.chi2s; % For readParameterAnnotation in fkt stringListChooser.m
    else
        ar.ple = [];
    end
else
    ar.ple = [];
end

save([savepath filesep, filename],'ar','-v7.3');
