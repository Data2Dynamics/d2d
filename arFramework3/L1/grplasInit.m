% GRPLASINIT Necessary initialization for Group Lasso scanning routine
%   Sets the relevant properties of global ar
% 
% grplasInit(jks, [means], [lbs], [ubs], linv, grouping, [weights], [thresh], [refit] )
% 
%   jks     parameter indices that shall be penalized
%   means   mean values of penalized parameters 
%           [0]
%   lbs     lower bounds 
%           [-5]
%   ubs     upper bounds 
%           [+5]
%   linv    inverse penalization strength vector to be scanned over
%   grouping determines the collection of parameters:
%           - integer vector in the size of jks: manual setting of groups
%           - 'single': each parameter has its own group -> REFER TO L1
%           instead
%           - 'altX' : with X integer, sets X fold alternate grouping. For jks =
%           [1 2 3 4] and grouping = 'alt2', group 1 will contain parameter
%           1 and 3, group 2 number 2 and 4
%           - 'groupX' : with X integer, clusters X subsequent parameters
%           together. For jks = [1 2 3 4] and grouping = 'group2', group 1
%           contains #1 and #2, group 2 contains number #3 and #4.
%   weights [1]
%           determines the weights applied to each group:
%           - 1x1 float: same weight for all groups
%           - length(jks) float array: sets a specific weight for each
%           parameter
%           - length(jks) / X float array with X referring to grouping
%           integer as above, then one weight is set per group
%           - [jks jks] size matrix: manual matrix
%           - [GL GL] size matrix with GL length of some group will assign
%           this matrix to each group of size GL. Other groups will have
%           the identity.
%   thresh  [1e-6]
%           determines the threshold epsilon below which a parameter is
%           considered to be equal to zero
%   refit   [true]
%           This flag indicates whether the model is fitted prior to
%           the analysis.
%           Refit is required if there are priors for the scanned
%           parameters since these priors are overwritten (temporarily)
%           for the analysis by L1-priors
% 
% See also l1Init, grplasScan


function grplasInit( jks, means, lbs, ubs, linv, grouping, weights, thresh, refit )
global ar

if(isempty(ar))
    error('please initialize by arInit')
end

ar.grplas = struct;

if(~exist('jks','var') || isempty(jks))
    arPrint
    fprintf('\n')
    error('please specify L1 parameters from above')
end

if(~exist('linv','var') || isempty(linv))
    error('please specify lambda scan range already in grplasInit')
else
    ar.linv = linv;
end

if(~exist('means','var') || isempty(means))
    means = zeros(size(jks));
end

if(~exist('lbs','var') || isempty(lbs))
    lbs = ones(size(jks))*(-5);
end

if(~exist('ubs','var') || isempty(ubs))
    ubs = ones(size(jks))*5;
end

if(~exist('weights','var') || isempty(weights))
    weights = 1;
end
        
if(~exist('thresh','var') || isempty(thresh))
    thresh = 1e-6;
end

if(~exist('refit','var') || isempty(refit))
    refit = true;
end

try
    if refit
        ar.qFit(jks) = 1;
        ar.type(jks) = 0;
        tmpar = ar;
        try
            arFit(true)
        catch
            ar = tmpar;
            arFitLHS(50)
        end
    end
    ar.grplas.lam0chi2s = arGetMerit('chi2')...
        ./ar.config.fiterrors_correction...
        +arGetMerit('chi2err');
end


ar.mean(jks) = means;
ar.lb(jks) = lbs;
ar.ub(jks) = ubs;

ar.type(jks) = 5;
ar.qFit(jks) = 1;
ar.std(jks) = Inf;

grplasSetGrouping(grouping,jks,weights);

ar.grplas.thresh = thresh;

end

