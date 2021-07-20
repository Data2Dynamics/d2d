%==========================================================================
% Matlab toolbox for structural identifiability and observability 
% analysis of nonlinear models
%--------------------------------------------------------------------------
% the toolbox based on:
% StrIkE-GOLDD (v3.0, last modified: 19/10/2020)
% https://github.com/afvillaverde/strike-goldd
%--------------------------------------------------------------------------


function arSIA(varargin)

global ar

% auto initialization
if ~exist('ar','var') || ~isfield(ar,'ia')
    fprintf(2, '\nThe analysis should be initialized by arSIAInit. It will be initialized with default settings \n');
    arSIAInit
end

% save the current path
currentFolder = pwd;

% add paths to the matlab search path
addpath(ar.ia.addedpaths{:})

% run the toolbox
STRIKE_GOLDD()

% auto reparametrization
if ar.ia.opts.autoRepar == 1
    fprintf('\n Auto reparametrization is starting\n');
    AutoRepar
end

% remove paths from the matlab search path
rmpath(ar.ia.addedpaths{:})

% redirect to the original path
cd(currentFolder)


end

