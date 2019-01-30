% ps = arRandomPars([n, randomseed])
% generate random parameter samples
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
% 
%   Only fitted parameters are sampled, i.e. only parameters with
%       a) qFit==1
%       b) qError==1 if ar.config.fiterrors==1
%   
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)

function ps = arRandomPars(n, randomseed)

global ar

if(~exist('n','var'))
    n = 10;
end

% set random seed
if(exist('rng','file')~=0)
    if(exist('randomseed','var') && ~isempty(randomseed))
        ar.lhs_seed = randomseed;
        rng(randomseed);
    else
        rng('shuffle');
        rngsettings = rng;
        ar.lhs_seed = rngsettings.Seed;
    end
end

% default matrix
ps = ones(n,1) * ar.p(:)';

lb = ar.lb;
ub = ar.ub;

% Use different bounds for initial guesses?
if isfield( ar, 'fitlb' )
    if ( numel( ar.fitlb ) ~= numel( ar.lb ) )
        error( 'ar.fitlb not the same size as ar.lb' );
    end
    lb = ar.fitlb;
    disp( 'Using alternate lower bounds' );
end
if isfield( ar, 'fitub' )
    if ( numel( ar.fitub ) ~= numel( ar.ub ) )
        error( 'ar.fitub not the same size as ar.ub' );
    end
    ub = ar.fitub;
    disp( 'Using alternate upper bounds' );
end

errorFitting = ( ar.config.fiterrors == 1) || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 );
if(isfield(ar.config, 'useLHS') && ar.config.useLHS==1) % LHS samples
    if errorFitting
        q_select = ar.qFit==1;
    else
        q_select = ar.qFit==1 & ar.qError~=1;        
    end
    
    psrand = lhsdesign(n,sum(q_select));
    psrand = psrand .* (ones(n,1)*(ub(q_select) - lb(q_select)));
    psrand = psrand + (ones(n,1)*lb(q_select));
    ps(:,q_select) = psrand;
elseif(isfield(ar.config, 'useLHS') && ar.config.useLHS==2) % random samples without LHS, prior considered if available
    for jp=1:length(ar.p)
        if( ar.qFit(jp)==1 && ~( ~errorFitting && ar.qError(jp)==1 ) )  % Error parameters should not be altered if we are not fitting errors
            if(ar.type(jp)==0 || ar.type(jp)==2) % uniform prior or uniform with normal bounds
                ps(:,jp) = lb(jp) + (ub(jp) - lb(jp)) * rand(n,1);
            elseif(ar.type(jp)==1 || ar.type(jp)==3) % normal prior or L1
                psrand = ar.mean(jp) + ar.std(jp) * randn(n,1);
                psrand(psrand>ub(jp)) = ub(jp);
                psrand(psrand<lb(jp)) = lb(jp);
                ps(:,jp) = psrand;
            else
                error('unsupported prior type');
            end
        end
    end
else % uniformly distributed, i.e. rand within the range [ar.lb, ar.ub]
    if errorFitting
        q_select = ar.qFit==1;
    else
        q_select = ar.qFit==1 & ar.qError~=1;
    end
    
    psrand = rand(n,sum(q_select));
    tmpublb = (ub(q_select) - lb(q_select));
    tmplb = lb(q_select);
    psrand = psrand .* (ones(n,1) * tmpublb(:)');
    psrand = psrand + (ones(n,1) * tmplb(:)');
    ps(:,q_select) = psrand;
end

