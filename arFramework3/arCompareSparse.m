%  Calculates the time required for arChi2 for
%     ar.config.useSparseJac=0; 
%   and
%     ar.config.useSparseJac=1;
% 
%   niter       number of calls of arChi2
%   ptrial_opt  0   use always ar.p
%               1   rand between ar.lb and ar.ub
%               2   randn around -1 with SD=1.
% 
% Useful Examples:
% [tnormal,tsparse] = arCompareSparse;
% [tnormal,tsparse] = arCompareSparse(20,1);
% [tnormal,tsparse] = arCompareSparse(20,2);


function [tnormal,tsparse] = arCompareSparse(niter,ptrial_opt)
if ~exist('niter','var') || isempty(niter)
    niter = 100;
end
if ~exist('ptrial_opt','var') || isempty(ptrial_opt)
    ptrial_opt = 1;
end
global ar

pIn = ar.p;

try
    
    tsparse = 0;
    tnormal = 0;
    
    for i=1:niter
        ar.config.useSparseJac=1;
        switch ptrial_opt
            case 0
                ptrial = ar.p;
            case 1  % rand between lb and ub
                ptrial = rand(size(ar.p)).*(ar.ub-ar.lb)+ar.lb;
            case 2  % randn around -1, std=1
                ptrial = randn(size(ar.p))-1;
        end
        
        tic
        arChi2(true,ptrial)
        tsparse = tsparse+toc;
        
        ar.config.useSparseJac=0;
        tic
        arChi2(true,ptrial)
        tnormal = tnormal+toc;
    end
catch ERR
    ar.p = pIn;
    rethrow(ERR)
end
ar.p = pIn;

tnormal = tnormal/niter;
tsparse = tsparse/niter;

fprintf('\nAverage tsparse = %.3f, tnormal = %.3f, tsparse/tnormal = %.3f\n\n',tsparse, tnormal, tsparse/tnormal);
