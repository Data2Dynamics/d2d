% Initialization of L1 scan
% jks    relative parameters to be investigated by L1 regularization

function l1Init(jks, means, lbs, ubs, linv, thresh, type, range, user_estim,refit)

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

if(~exist('type','var') || isempty(type))
    type = 'L1';
elseif(~exist('range','var') || isempty(range))
    type = 'L1';
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



if (strcmpi(type,'L1'))
    ar.L1subtype(jks) = 1;
    % classical L1 penalization, no further requirements
elseif (strcmpi(type,'adaptive') || strcmpi(type,'al'))
    
    
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
    
    
    
elseif (strcmpi(type,'Lq') || strcmpi(type,'Lnu'))
    
    
    ar.L1subtype(jks) = 3;
    ar.nu = zeros(size(linv)) + range;
    ar.expo(jks) = ar.nu(1);
    % range must be a scalar or a linv-size vector
    
elseif (strcmpi(type,'elasticnet') || strcmpi(type,'en'))
    
    ar.L1subtype(jks) = 4;
    ar.alpharange = zeros(size(linv)) + range;
    ar.alpha(jks) = ar.alpharange(1);
    
    
end








ar.type = oldtype;
ar.type(jks) = 3;
ar.qFit(jks) = 1;
ar.std(jks) = Inf;



% try
%     arPrint(jks)
% end
    
ar.L1ps = [];
ar.L1chi2s = [];
ar.L1chi2fits = [];

ar.L1ps_unpen = [];
ar.L1chi2s_unpen = [];
ar.L1thresh = thresh;