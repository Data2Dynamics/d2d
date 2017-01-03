% ok = arCheckSimu(whichp, varargin)
% 
% This function checks whether arSimu is feasible for the current
% parameters. If not, the outcome is checked if parameter
% components are set to -1. 
% 
%   whichp  specifies which parameters should be tested, i.e. set to -1
%           default: whichp = find(ar.qDynamic)
% 
%   varargin Arguments passed to arSimu
% 
%   ok      1: Integration is feasible for ar.p
%           [0,1,0,1,1,NaN,1]: zero indicates non-feasibility of integration
%           if parameter is set to -1
%           one indicates that feasibitly
%           NaN indicates that the parameter was not checked
% 
% The function is desired to find out which parameters are responsible for
% convergence failure during ODE integration.


function ok = arCheckSimu(whichp, varargin)
global ar

if(~exist('whichp','var') | isempty(whichp))
    whichp = find(ar.qDynamic);
elseif(size(whichp,1)>1)
    whichp = whichp'; %should not be a column
end

if(ischar(whichp))
    whichp = find(~cellfun(@isempty,regexp(ar.pLabel,whichp)));
elseif islogical(whichp)
    whichp = find(whichp);
elseif iscell(whichp)
    js2 = NaN(size(whichp));
    for i=1:length(whichp)
        tmp = strmatch(whichp{i},ar.pLabel,'exact');
        if length(tmp)==1
            js2(i) = tmp;
        end
    end
    whichp = js2;
elseif(sum(isnan(whichp))>0 || sum(isinf(whichp))>0 || min(whichp)<1 || max(whichp-round(whichp))>eps)
    whichp
    warning('arCheckSimu.m: argument whichp is not plausible.')
end

pIn = ar.p + 0.0; % +0.0 because of MATLAB bug in some releases 
inok = 0;
try 
    arSimu(varargin{:})
    inok = 1;
end


if(~inok)
    disp('Failure for ar.p with the following message:');
    disp(lasterr)
    
    ok = NaN(size(ar.p));
    ok(whichp) = 0;
    for i=1:length(whichp)
        try
            ar.p(whichp(i)) = -1;
            fprintf(['Testing %12s\t = -1: '],ar.pLabel{whichp(i)});
            arSimu(varargin{:})
            fprintf(' -> Integration is now feasible! This parameter seems to cause integration problems. \n')

            ar.p = pIn + 0.0; % +0.0 because of MATLAB bug in some releases
            ok(whichp(i))=1;
        catch
            fprintf(' Integration not yet feasible.\n')
%             disp(lasterr)
            ar.p = pIn + 0.0; % +0.0 because of MATLAB bug in some releases
        end
    end

    if(sum(ok==1)>0)
        fprintf('\nThe following parameters are responsible for the integration errors:\n')
        fprintf('%s ',ar.pLabel{find(ok==1)});
        fprintf('\n\n')
    end
else
    if nargout==0
        fprintf('Integration with the initial parameters (ar.p), feasible => Checking the parameters for causing intergration problems omitted.');
    end
    ok = true;
end
