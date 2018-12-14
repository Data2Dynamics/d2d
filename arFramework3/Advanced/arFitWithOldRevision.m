% chi2 = arFitWithOldRevision(sha,varargin)
%
% This script refits using the original D2D version. This only works, if
% ar.info.revision contains the sha (which only works if D2D is used in
% combination with git).
%
%   sha         ID of the different commits/revisions on github
%               (eg 7182eff7eb6f723da7b45da5470a37fd13c9e1de)
%               see commit pages at github
%               the commit sha is stored in ar.info.revision
%               [ar.info.revision]
%   varargin    fitting options passed to arFit(varargin{:})
%
%   chi2        ar.chi2fit
% 
% Attention: This function does not recompile using the old revision. It
% uses/requires existing mex file as specified via ar.fkt
% 
% If you want to compare Revision including recompilation, use
% arEvalWithOldRevision.m e.g. via
% arEvalWithOldRevision(@arFitLHS,'25b298f0957ddafbccf04d805593ef02355bb3dc',true)
%
% Example:
%   arFitWithOldRevision('c949d75ad3a00ee92248cfaffb6060d5c01889e8')


function chi2 = arFitWithOldRevision(sha,varargin)
global ar

if ~exist('sha','var') || isempty(sha)
    if isfield(ar.info,'revision')
        sha = ar.info.revision;
    else
        disp('Revision unknown, may because you are not using git.');
        return;
    end
end


revision_path = arGetOldRevision(sha);

%% set path to old revision:
addpath([revision_path,filesep,'d2d-',sha,filesep,'arFramework3']);
ar_path = fileparts(which('arInit.m'));
fprintf('ar_path is now %s \n.',ar_path);

tmp_paths = genpath(ar_path);
addpath(tmp_paths);


%% Do fitting
try
    arFit(varargin{:})
    if exist([revision_path,filesep,'d2d-',sha,filesep,'arFramework3',filesep,'arGetMerit'],'file')
        chi2 = arGetMerit('chi2fit');
    else % older versions
        chi2 = ar.chi2fit;
    end

    rmpath(tmp_paths);
    fprintf('Fit with revision %s finished.\n',sha);
catch ERR
    ar.fkt
    fprintf('Fitting with an old revision failed.\n');        
    rmpath(tmp_paths);
    
    rethrow(ERR)
end
