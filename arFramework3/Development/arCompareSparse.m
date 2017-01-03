%  Calculates the time required for arChi2 for
%     ar.config.useSparseJac=0; 
%   and
%     ar.config.useSparseJac=1;
% 
%   niter       number of calls of arChi2
%   ptrial_opt  0   use always ar.p
%               1   rand between ar.lb and ar.ub
%               2   randn around -1 with SD=1.
%               3   parameters provided as argument ptrial_in
% 
%   fun_opt  1   arChi2
%            2   arFit
% 
%   ptrial_in   each row is a parametre vector used for assemssment if
%               ptrial_opt==3
% 
% Useful Examples:
% [tnormal,tsparse] = arCompareSparse;
% [tnormal,tsparse] = arCompareSparse(20,1);     % p drawn with rand, arChi2 
% [tnormal,tsparse] = arCompareSparse(20,2);     % p drawn with randn, arChi2
% [tnormal,tsparse] = arCompareSparse(20,2,2);   % p drawn with randn, arFit 
% [t0,t1]=arCompareSparse2(10,3,1,ar.fit.p_hist) % following a path from the previous fit, arChi2

function [tnormal,tsparse] = arCompareSparse(niter,ptrial_opt,fun_opt,ptrial_in)
if ~exist('niter','var') || isempty(niter)
    niter = 100;
end
if ~exist('ptrial_opt','var') || isempty(ptrial_opt)
    ptrial_opt = 1;
end
if ~exist('ptrial_in','var') || isempty(ptrial_in)
    ptrial_in = [];
end
if ~exist('fun_opt','var') || isempty(fun_opt)
    fun_opt = 1; % arChi2
end

global ar

pIn = ar.p;

if(ptrial_opt==3)
    if isempty(ptrial_opt)
        error('ptrial_opt=3 requires a matrix of pameters ptrial_in as third argument');
    else
        niter = sum(sum(isnan(ptrial_in),2)==0,1);
    end
end


try    
    tsparse = 0;
    tnormal = 0;
    
    ptrial = NaN(niter,sum(ar.qFit==1));
    for i=1:niter        
        switch ptrial_opt
            case 0
                ptrial(i,:) = ar.p(ar.qFit==1);
            case 1  % rand between lb and ub
                ptrial(i,:) = rand(size(ar.p(ar.qFit==1))).*(ar.ub(ar.qFit==1)-ar.lb(ar.qFit==1))+ar.lb(ar.qFit==1);
            case 2  % randn around -1, std=1
                ptrial(i,:) = randn(size(ar.p(ar.qFit==1)))-1;
            case 3
                ptrial(i,:) = ptrial_in(i,ar.qFit==1);
        end
    end
    
        
    tsparse = NaN(niter,1);
    tnormal = NaN(niter,1);
    for i=1:niter
        ar.config.useSparseJac=0;
        switch fun_opt
            case 1
                tic
                arChi2(true,ptrial(i,:))
            case 2
                ar.p(ar.qFit==1) = ptrial(i,ar.qFit==1);
                tic
                arFit
        end
        tnormal(i)=toc;
%     end
%     
%     for i=1:niter
        ar.config.useSparseJac=1;        
        switch fun_opt
            case 1
                tic
                arChi2(true,ptrial(i,:))
            case 2
                ar.p(ar.qFit==1) = ptrial(i,ar.qFit==1);
                tic
                arFit
        end
        tsparse(i) = toc;

        if fun_opt==2
            fprintf('tnormal = %f, tsparse=%f\n\n',tnormal(i),tsparse(i));
        end
    end
catch ERR
    ar.p = pIn;
    rethrow(ERR)
end
ar.p = pIn;

fprintf('\nmean(tsparse) = %f, mean(tnormal) = %f, mean(tsparse)/mean(tnormal) = %.3f\n\n',mean(tsparse), mean(tnormal), mean(tsparse)/mean(tnormal));

