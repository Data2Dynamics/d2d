% l1Init(jks, [means], [lbs], [ubs], linv, [thresh], [regtype], (range), user_estim, refit)
%
% Initialization of L1 scans
% 
% jks           parameter indices of the investigated parameters by L1 regularization
%               Typically these parameters denot fold-factors.
% means         Target values of the L1-priors, i.e. prior means
% lbs           Lower bounds                [-5*ones(size(jks))]
% ubs           Upper bounds                [5*ones(size(jks))]
% linv          The range for the inverse penalization parameter lambda
%               [ flip(logspace(-4,0,40)) ]
% thresh        [1e-6]
%               The threshold for deciding whether a parmeter conicides
%               with the L1 target value ar.mean. This threshold is stored
%               in ar.L1thresh. 
% regtype       ['L1']      standard Lasso
%               'adaptive'  adaptive Lasso
%               'Lq'
%               'elasticnet'
% range         []
%               Has to be specified if regtype is not 'L1'
%               range is added to 
%               - ar.gamma for adaptive L1
%               - ar.nu for Lq
%               - ar.alpha for elasticnet
% user_estim    [ar.p(jks)]
%               user-defined estimates, only used for adaptive lasso for
%               calculating ar.lnuweights
% refit         [false]
%               This flag indicates whether the model is fitted prior to
%               the analysis.
%               Refit is required if there are priors for the scanned
%               parameters since these priors are overwritten (temporarily)
%               for the analysis by L1-priors


function l1Init(jks, means, lbs, ubs, linv, thresh, regtype, range, user_estim, refit)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    arPrint
    fprintf('\n')
    error('please specify L1 parameters from above')
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

if(~exist('linv','var') || isempty(linv))
    warning('please specify lambda scan range in l1Init')
    linv = logspace(-4,0,40);
    linv = linv(end:-1:1);
end

ar.mean(jks) = means;
ar.lb(jks) = lbs;
ar.ub(jks) = ubs;
ar.linv = linv;



if(~exist('thresh','var') || isempty(thresh))
    thresh = 1e-6;
end

if(~exist('regtype','var') || isempty(regtype))
    regtype = 'L1';
end
if(~exist('range','var') || isempty(range))
    if strcmpi(regtype,'l1')~=1
        error('Input argument range has to be specified for non-L1 penalization');
    end
end

if(~exist('refit','var') || isempty (refit))
    refit = false;
end

% Unpenalized chi2 values for comparison against the penalized
% ones later on in the SELECTOPT procedure
% In default setting, the OLS estimation is used for weighting the L1

oldtype = ar.type;

try
    if refit
        ar.qFit(jks) = 1;
        ar.type(jks) = 0;
        arFit(true)
    end
    ar.L1lam0chi2s = arGetMerit('chi2')./ar.config.fiterrors_correction...
        +arGetMerit('chi2err');
    estim = ar.p(jks);
    estim(estim == 0) = 1e-15;
end



if (strcmpi(regtype,'L1'))
    ar.L1subtype(jks) = 1;
    % classical L1 penalization, no further requirements

elseif (strcmpi(regtype,'adaptive') || strcmpi(regtype,'al'))    
    
    ar.L1subtype(jks) = 2;
    ar.gamma = zeros(size(linv)) + range;
    % range must be a scalar or a linv-size vector
    if (~exist('user_estim','var') || isempty(user_estim) || strcmpi(user_estim, 'OLS'))
        ar.estim(jks) = estim;
        % default is OLS estimation calculated above
    else
        ar.estim(jks) = user_estim;
        % user provided estimations
    end
    ar.lnuweights(jks) = 1./ar.estim(1).^ar.gamma(1);
        
elseif (strcmpi(regtype,'Lq') || strcmpi(regtype,'Lnu'))
        
    ar.L1subtype(jks) = 3;
    ar.nu = zeros(size(linv)) + range;
    ar.expo(jks) = ar.nu(1);
    % range must be a scalar or a linv-size vector
    
elseif (strcmpi(regtype,'elasticnet') || strcmpi(regtype,'en'))
    
    ar.L1subtype(jks) = 4;
    ar.alpharange = zeros(size(linv)) + range;
    ar.alpha(jks) = ar.alpharange(1);
    
else
    error('Regularization type %s is not implemented.',regtype);
end

if sum(ar.type(jks)~=0)>0
    pnames = sprintf('%s, ',ar.pLabel(jks));
    warning('ar.type and ar.std is overwritten for the scanned parameters: %s',pnames(1:end-2));
end

if ~isfield(ar, 'L1DiffPen_activate') || isempty(ar.L1DiffPen_activate)
    ar.L1DiffPen_activate = 0;
end

ar.type(jks) = 3;
ar.qFit(jks) = 1;
ar.std(jks) = Inf;
    
ar.L1ps = [];
ar.L1chi2s = [];
ar.L1chi2fits = [];

ar.L1ps_unpen = [];
ar.L1chi2s_unpen = [];
ar.L1thresh = thresh;