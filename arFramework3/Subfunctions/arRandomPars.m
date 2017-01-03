% generate random parameter samples
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
% 
%   Only fitted parameters are sampled, i.e. only parameters with
%       a) qFit==1
%       b) qError==1 if ar.config.fiterrors==1
%   
%
% ps = arRandomPars(n, randomseed)
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
ps = ones(n,1) * ar.p;

if(isfield(ar.config, 'useLHS') && ar.config.useLHS==1) % LHS samples
    if ar.config.fiterrors==1
        q_select = ar.qFit==1;
    else
        q_select = ar.qFit==1  & ar.qError~=1;
    end
    
    psrand = lhsdesign(n,sum(q_select));
    psrand = psrand .* (ones(n,1)*(ar.ub(q_select) - ar.lb(q_select)));
    psrand = psrand + (ones(n,1)*ar.lb(q_select));
    ps(:,q_select) = psrand;
elseif(isfield(ar.config, 'useLHS') && ar.config.useLHS==2) % random samples without LHS, prior considered if available
    for jp=1:length(ar.p)
        if(ar.qFit(jp)==1  && (ar.config.fiterrors~=0 || ar.qError(jp)~=1 ) )  % Error parameters should not be altered if fiterrors==0
            if(ar.type(jp)==0 || ar.type(jp)==2) % uniform prior or uniform with normal bounds
                ps(:,jp) = ar.lb(jp) + (ar.ub(jp) - ar.lb(jp)) * rand(n,1);
            elseif(ar.type(jp)==1 || ar.type(jp)==3) % normal prior or L1
                psrand = ar.mean(jp) + ar.std(jp) * randn(n,1);
                psrand(psrand>ar.ub(jp)) = ar.ub(jp);
                psrand(psrand<ar.lb(jp)) = ar.lb(jp);
                ps(:,jp) = psrand;
            else
                error('unsupported prior type');
            end
        end
    end
else % uniformly distributed, i.e. rand within the range [ar.lb, ar.ub]
    if ar.config.fiterrors==1
        q_select = ar.qFit==1;
    else
        q_select = ar.qFit==1  & ar.qError~=1;
    end
    
    psrand = rand(n,sum(q_select));
    psrand = psrand .* (ones(n,1)*(ar.ub(q_select) - ar.lb(q_select)));
    psrand = psrand + (ones(n,1)*ar.lb(q_select));
    ps(:,q_select) = psrand;
end
