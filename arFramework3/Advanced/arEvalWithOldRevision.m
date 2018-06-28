%   This script evaluatees a specific function using a defined D2D version.
%   The D2D revision is specified via the commit's shas (which only works if
%   D2D is used in combination with git).
%
%   dorecompile     Default: false
%
%                   With this option, the function does not recompile using
%                   the old revision. It uses/requires existing mex file as
%                   specified via ar.fkt
%
%                   true    if Setup.m is available, it is called.
%                   Otherwise, arRecompile is tried.
%
% Example:
% result = arEvalWithOldRevision; % Default shas without recompilation
% result = arEvalWithOldRevision([], [], true);  % with revision dependent compilation


function result = arEvalWithOldRevision(fun, shas, dorecompile, varargin)
global ar
if ~exist('fun','var') || isempty(fun)
    fun = @fitLHS; % This function is defined below (as template/default).
end

if ~exist('shas','var') || isempty(shas)
    shas = {'0bfcb5facea85773a262fe417a213e93b6681dfe',...  % 27.6.18
        'c949d75ad3a00ee92248cfaffb6060d5c01889e8',... % 31.5.18, before statistics toolbox was removed
        '8e5009ccbc51f5b6fce6db89c06ae9f60d738454',... % 6.3.18  new field for saving the arSimuCalc mex file
        '0e130ed7c9aa87902086b49ecf158ae7265c409b',... % 13.2.18,  Added relative tolerance to equilibration as well. Equilibration is c…
        '7182eff7eb6f723da7b45da5470a37fd13c9e1de',... % 20.11.2017  Included more descriptive message when equilibration fails
        'a60a059ee0d2672b7dc063ab9578b6e9f3575cce',... % BAD Becker, 9.10.2017  Fixes for custom functions with full expressions in their arguments.
        'b74a18a66205ef317579944ae3dcddcfba28541c',... % BAD Becker, 22.9.2017  Fixes rootfinding
        '25b298f0957ddafbccf04d805593ef02355bb3dc',... % 26.7.17 Big update: when setting ar.config.sensitivitySubset to 1, D2D now on…
        'd564b10acd435c16ba52ff635f10f244974f2e08',... % 15.5.17:  Added option to make it default to only save parameters with arSave
        '0d9217318ba24ce725e1a013a99aea2d5d6ad23e',... % 31.1.17:  Removed unused variables arSimuCalc.c, AFTER OBERSTDORF
        'ae24265ca6852c68dc65a50012574290cbd11b15',... % 21.12.16: Before Oberstdorf
        '5c3fc2251709ac958c21d7c2ed33c1ba2275a337',... % 30.12.2015:     added proper wait bar with time forecast to cluster functions
        };
elseif ischar(shas)
    shas = {shas};
end
if ~exist('dorecompile','var') || isempty(dorecompile)
    dorecompile = false;
end

rev0 = ar.info.revision;
fkt0 = ar.fkt;
    
result = cell(size(shas));
if dorecompile
    if exist([fkt0,'.',mexext],'file')
        % backup and remove ar.fkt (otherwise compilation might be skipped)
        my_munlock([fkt0,'.',mexext]) % to prevent blocking the function because "in use"
        movefile([fkt0,'.',mexext],[fkt0,'_backup.',mexext]);
    end
end

for s=1:length(shas)
    arRemoveOldRevisionPaths        
    
    revision_path = arGetOldRevision(shas{s});
    
    %% set path to old revision:
    addpath([revision_path,filesep,'d2d-',shas{s},filesep,'arFramework3']);
    ar_path = fileparts(which('arInit.m'));
    fprintf('ar_path is now %s \n.',ar_path);
    
    tmp_paths = genpath(ar_path);
    addpath(tmp_paths);
    
    %% recompilation?
    if dorecompile 
        if exist([ar.fkt,'_revision_',shas{s},'.',mexext],'file') % usually the combi of ar.fkt and shas ist not yet available before Setup
            my_munlock([ar.fkt,'_revision_',shas{s},'.',mexext]);
            my_munlock([ar.fkt,'.',mexext])
            copyfile([ar.fkt,'_revision_',shas{s},'.',mexext],[ar.fkt,'.',mexext]);
        elseif exist('Setup.m','file')
            disp('Setup will be executed ...')
            Setup;
            close all % setup often raise figures
        elseif exist('arRecompile','file')
            disp('arRecompile will be executed ...')
            arRecompile;
        else
            setupfile = input('Please specify Setup file: ','s');
            eval(setupfile);
            close all % setup often raise figures
        end
        
        if exist([ar.fkt,'.',mexext],'file')
            my_munlock([ar.fkt,'.',mexext])
            my_munlock([ar.fkt,'_revision_',shas{s},'.',mexext])
            copyfile([ar.fkt,'.',mexext],[ar.fkt,'_revision_',shas{s},'.',mexext]);  % make a copy for all revisions
        else
            error('%s not found. Check whether this function is (unintendedly) somewhere else in the Matlab path',[ar.fkt,'.',mexext]);
        end
    end
    
    
    %% evaluating the function
    try
        arCheckCache( true ) % just to be sure        
        ar.info.revision = shas{s};
        
        result{s} = feval(fun,varargin{:});        

        rmpath(tmp_paths);
        arRemoveOldRevisionPaths % to be sure
        ar.info.revision = rev0;
        if dorecompile % remove mex-file because you cannot see the revision (a copy incl. revision ID is made above)
            delete([ar.fkt,'.',mexext]);
        end
        ar.fkt = fkt0;

            
        fprintf('Fit with revision %s finished.\n',shas{s});
    
    catch ERR
        ar.fkt
        fprintf('Fit sequence %i, fit number %i could not be refitted. Maybe the mex file does not run on this system.\n');

        %%%%%%%%%%%%%%
        % reset revision:
        rmpath(tmp_paths);
        arRemoveOldRevisionPaths % to be sure
                
        ar.info.revision = rev0;
        ar.fkt = fkt0;
        if dorecompile
            if exist([fkt0,'_backup.',mexext],'var')
                movefile([fkt0,'_backup.',mexext],[fkt0,'.',mexext]);
            end
        end
        %%%%%%%%%%%%%%
        rethrow(ERR)
    end                
end

%%%%%%%%%%%%%%
% reset revision:
ar.info.revision = rev0;
ar.fkt = fkt0;
if dorecompile
    if exist([fkt0,'_backup.',mexext],'var')
        movefile([fkt0,'_backup.',mexext],[fkt0,'.',mexext]);
    end
end
%%%%%%%%%%%%%%

end

function res = fitLHS
global ar

load ps_start3

arFits(ps_start);

res = struct;
res.revision = ar.info.revision;
res.fkt = ar.fkt;
res.ar = arDeepCopy(ar);
res.ps = ar.ps;
res.chi2s = ar.chi2s;
res.ps_start = ar.ps_start;
end

function my_munlock(file)
try
    munlock(file)
end
end
