% arCompareFitsOverRevisions
% 
% arCompareFitsOverRevisions(shas,p0)
% 
%   shas    IDs specifying the commit
% 
%   p0      initial guess for the parameters
% 
% [pfit,chi2s,chi2new] = arCompareFitsOverRevisions
% 
%   pfit    fitted parameters
%   
%   chi2s   chi2fit according to the different revisions
% 
%   chi2new arGetMerit using the current revision (given by global ar)
% 
% 
% This function can be used to compare fit results over several D2D
% version. 
% Note that for the first time, this function is applied, the commits are
% downloaded via arGetOldRevision which takes time.
% 
% 
% Example:
% [pfit,chi2s] = arCompareFitsOverRevisions

function [pfit,chi2s,chi2new] = arCompareFitsOverRevisions(shas,p0)

if ~exist('shas','var') || isempty(shas)
    shas = {'0bfcb5facea85773a262fe417a213e93b6681dfe',...  % 27.6.18
        'c949d75ad3a00ee92248cfaffb6060d5c01889e8',... % 31.5.18, before statistics toolbox was removed
        '8e5009ccbc51f5b6fce6db89c06ae9f60d738454',... % 6.3.18  new field for saving the arSimuCalc mex file
        '0e130ed7c9aa87902086b49ecf158ae7265c409b',... % 13.2.18,  Added relative tolerance to equilibration as well. Equilibration is c…
        '7182eff7eb6f723da7b45da5470a37fd13c9e1de',... % 20.11.2017  Included more descriptive message when equilibration fails
        'a60a059ee0d2672b7dc063ab9578b6e9f3575cce',... % 9.10.2017  Fixes for custom functions with full expressions in their arguments.
        'b74a18a66205ef317579944ae3dcddcfba28541c',... % 22.9.2017  Fixes rootfinding
        '25b298f0957ddafbccf04d805593ef02355bb3dc',... % 26.7.17 Big update: when setting ar.config.sensitivitySubset to 1, D2D now on…
        'd564b10acd435c16ba52ff635f10f244974f2e08',... % 15.5.17:  Added option to make it default to only save parameters with arSave
        '0d9217318ba24ce725e1a013a99aea2d5d6ad23e',... % 31.1.17:  Removed unused variables arSimuCalc.c, AFTER OBERSTDORF
        'ae24265ca6852c68dc65a50012574290cbd11b15',... % 21.12.16: Before Oberstdorf
        '5c3fc2251709ac958c21d7c2ed33c1ba2275a337',... % 30.12.2015:     added proper wait bar with time forecast to cluster functions
        };  
elseif ~iscell(shas) % if only a single sha is provided 
    shas = {shas};
end

global ar

if ~exist('p0','var') || isempty(p0)
    p0 = ones(size(ar.p));
end

pfit  = cell(size(shas));
chi2s = NaN(size(shas));

try 
    delete('arFitWithOldRevision.log') % enforce deleting old stuff
end
diary('arFitWithOldRevision.log');
for i=1:length(shas)
    ar.p = p0;
    
    diary off % enforce writing, because matlab might crash for outdated c-code which will be called for old revisions
    diary('arFitWithOldRevision.log');
    
    fprintf('Fit with revision %s now ...\n',shas{i});
    chi2s(i) = arFitWithOldRevision(shas{i});
    pfit{i} = ar.p;    
    chi2s
end

chi2new = NaN(size(chi2s)); % Objective function according to the current (given by the input) code.
for i=1:length(pfit)
    ar.p = pfit{i};
    arCalcMerit;
    chi2new(i) = arGetMerit(true);
end


