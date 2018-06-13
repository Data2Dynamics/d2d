% MC^3 Sampler (contact: franz-georg.wieland@mars.uni-freiburg.de)
% MCMC sampler with multiple sampling methods and Parallel Tempering
% based on previous code by Joep Vanlier
%               and code by Benjamin Ballnus in PESTO
% inspired by code by M. Girolami 
%
% function arMCMC(nruns, nburnin, method, exchange)
%
%   nruns 
%               Number of runs
%   nburnin
%               Number of Burn-In iterations to be discarded (only for
%              adaptive MCMC)
%   method for proposal density:
%           1 = N(0,1)
%           2 = N(0,c) scaled
%           3 = Adaptive MCMC
%           4 = Fisher based
%           5 = Simplified mMALA SVD
%           6 = Simplified mMALA Cholesky
%   exchange
%               Thinning coefficient (100 means that only every 100th point
%               is used for posterior)
%               Increase to reduce autocorrelation between samples
%  
%  Normal Distributions
%   optimal acceptance rate ~23% (for multi-variate normal distributions),
%   see in:
%   Roberts, G.O.; Gelman, A.; Gilks, W.R. (1997).
%   "Weak convergence and optimal scaling of random walk Metropolis algorithms".
%   Ann. Appl. Probab. 7 (1): 110-120.
%
%  MMALA
%   Girolami, M. and Calderhead, B. (2011).
%   "Riemann manifold Langevin and Hamiltonian Monte Carlo methods"
%   JRSS B 73 (2): 123-214.
%   Optimal acceptance rate for MMALA ~57.4%. See
%   OPTIMAL SCALING AND DIFFUSION LIMITS FOR THE
%   LANGEVIN ALGORITHM IN HIGH DIMENSIONS
%   Natesh S. Pillai, Andrew M. Stuart and Alexandre H. Thiery
%   Ann. Appl. Probab. 22 (6): 2320-2356

function mcmc = arMC3(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Readout and Initialization of global function variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global ar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Readout ar.mc3 struct for specified options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manual Scaling of Covariance Matrix to account for higher dimensions

if isfield(ar, 'mc3')

    if isfield(ar.mc3, 'nruns')
        nruns = ar.mc3.nruns;
    else
        nruns = 10000;
    end
    
    if isfield(ar.mc3, 'nburnin')
        nburnin = ar.mc3.nburnin;
    else
        nburnin = 0;
    end
    
    if isfield(ar.mc3, 'method')
        method = ar.mc3.method;
    else
        method = 4;
    end
    
    if isfield(ar.mc3, 'exchange_method')
        exchange_method = ar.mc3.exchange_method;
    else
        exchange_method = 0;
    end       
    
    if isfield(ar.mc3, 'nthinning')
        nthinning = ar.mc3.nthinning;
    else
        nthinning = 1;
    end
    
    if isfield(ar.mc3, 'ManualScalingFactor')
        ManualScalingFactor = ar.mc3.ManualScalingFactor;
    else
        ManualScalingFactor = 1;
    end
    
    if isfield(ar.mc3, 'UseScaling')
        UseScaling = ar.mc3.UseScaling;
    else
        UseScaling = 1;
    end
    
    if isfield(ar.mc3, 'TemperatureExponent')
        TemperatureExponent = ar.mc3.TemperatureExponent;
    else
       TemperatureExponent = 1;
    end  

    if isfield(ar.mc3, 'RegularizationThreshold')
        RegularizationThreshold = ar.mc3.RegularizationThreshold;
    else
        RegularizationThreshold = 1e-8;
    end
    
    if isfield(ar.mc3, 'min_accept')
        min_accept = ar.mc3.min_accept;
    else
        min_accept =  0.4;
    end
    
    if isfield(ar.mc3, 'max_accept')
        max_accept = ar.mc3.max_accept;
    else
       max_accept = 0.7;
    end  
%     
%     if isfield(ar.mc3, 'refresh_rate')
%         refresh_rate = ar.mc3.refresh_rate;
%     else
%         refresh_rate = 500;
%     end  
% 
%     if isfield(ar.mc3, 'show')
%         show = ar.mc3.show;
%     else
%         show = 0;
%     end
    
    if isfield(ar.mc3, 'Cmax')
        Cmax = ar.mc3.Cmax;
    else
        Cmax = 1e8;
    end
    
    if isfield(ar.mc3, 'Cmin')
        Cmin = ar.mc3.Cmin;
    else
       Cmin = 1e-15;
    end      
    
    if isfield(ar.mc3, 'Cmod')
        Cmod = ar.mc3.Cmod;
    else
       Cmod = 1.01;
    end    
    
    if isfield(ar.mc3, 'NumberOfChains')
        NumberOfChains = ar.mc3.NumberOfChains;
    else
       NumberOfChains = 1;
    end    
    
    if isfield(ar.mc3, 'LengthOfAcceptanceTestChain')
        LengthOfAcceptanceTestChain = ar.mc3.LengthOfAcceptanceTestChain;
    else
        LengthOfAcceptanceTestChain          = 25;
    end
    
    if isfield(ar.mc3, 'DecayParameter')
        DecayParameter = ar.mc3.DecayParameter;
    else
        DecayParameter = 0.51;
    end    

end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % Set defaults
% [nruns, nburnin, method, exchange_method] = deal( ...
%     10000,...                % nruns         Number of runs
%     0,...                   % nburnin       Number of burn-in runs
%     5,...                   % method
%     0 ...                   % exchange      Choice of exchange function
%     );

% Read out given variables and overwrite struct varibales if given
if nargin > 3
	exchange_method = varargin{4};
end
if nargin > 2
	method = varargin{3};
end
if nargin > 1
	nburnin = varargin{2};
end
if nargin > 0
    nruns = varargin{1};    
end

% Counter for Number of Bound Violations
NumberOfBoundViolations    = 0;



if exchange_method == 0 && NumberOfChains >1
    NumberOfChains = 1;
    disp( 'Warning: Number of Chains is set to 1 because no exchange method was choosen!');
    disp( '         Use ar.mc3.exchange_method to choose a exchange method.');
end


beta = linspace(1,1/NumberOfChains,NumberOfChains).^TemperatureExponent;    % Scaling of Temperature regime based on 


paraReset            = repmat( ar.p, 1, 1, NumberOfChains );               % Parameter Reset vector with current parameter value for all chains
qFitGlobal           = ar.qFit==1;                                   % Logical variable - true for fitted parameters, false for fixed as specified in Global model input
ResidualType         = ar.res_type;                                  % Type of the residual: 1 (Data point residual),2 (Residuals of fitted errors, 3 (Prior), 4 (Random Effects)
locs = 1:sum((qFitGlobal));



% Define nwindow as 50 times number of fitted parameters
% Window will define the number of parameter samples that are used in
% adaptive MCMC for adjusting during Burn-In phase
nwindow = sum(qFitGlobal)*50;

% Creation of struct elements t it does not yet exist    
ar.mc3.nwindow = nwindow;    

% Initialize storing variables
parasHistory = nan(nwindow, sum(qFitGlobal), NumberOfChains);
parasHistory_index = 1;





ProbProposalGivenCurr = 1e-30.*ones(NumberOfChains,1);
ProbCurrGivenProposal = 1e-30.*ones(NumberOfChains,1);

% Cfactor = (2.38/sqrt(sum(qFit)))^2 / sum(qFit);
Cfactor = Cmin*Cmax * ones(NumberOfChains,1);

% bounds
lb = ar.lb(qFitGlobal);
ub = ar.ub(qFitGlobal);



% Check if Error Correction Fit is enables, if so disable it because this
% is not used in Bayesian Setting
if isfield( ar, 'useFitErrorCorrection' )
    fitErr = ar.useFitErrorCorrection;
    if ( fitErr )
        disp( 'Warning: Use Fit Error Correction was on. Disabled for MCMC since');
        disp( 'we need the actual likelihood function!' );
        ar.useFitErrorCorrection = 0;
    end
else
    fitErr = -1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize chains and iteratively updated variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para_curr = repmat( ar.p(qFitGlobal), NumberOfChains, 1 );

% % initial values for chains
% if (isfield(ar,'mcmc'))
%     if (isfield(ar.mcmc,'p'))
%         % Multiple initial parameter vector supplied
%         para_curr = ar.mcmc.p(:,qFitGlobal);
%     else
%         if ( NumberOfChains > 1 )
%             disp( 'Warning: All chains starting from same parameter vector' );
%         end
%         para_curr = repmat( ar.p(qFitGlobal), NumberOfChains, 1 );
%     end
% else
%     para_curr = repmat( ar.p(qFitGlobal), NumberOfChains, 1 );
% end    


jindexoffset = 0;
ar.ps = nan(nruns, length(ar.p), NumberOfChains);
ar.ps_trial = nan(nruns, length(ar.p), NumberOfChains);
ar.Q_trial = nan(nruns,NumberOfChains);
ar.Q_curr = nan(nruns,NumberOfChains);
ar.chi2s = nan(nruns,NumberOfChains);
ar.chi2s_trial = nan(nruns,NumberOfChains);
ar.acceptance = nan(nruns,NumberOfChains);
ar.exchange = nan(nruns,NumberOfChains);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up chosen sampling functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Adaptive MCMC has a minimal burn in requirement
% nburnin = max( [ nburnin, (method==4)*nwindow ] );

nburnin = max( [ nburnin, (method==3)*nwindow, (method==2)*nwindow  ] );

% methods
switch method
    case 1
        mName = 'Unscaled MVN MCMC';
        fkt = @mcmc_norm_mvnrnd;
        adjust_scaling = false;
        use_sensis = false;
        burnin = false;
    case 2
        mName = 'Scaled MVN MCMC';
        fkt = @mcmc_scaled_mvnrnd;
        adjust_scaling = true;
        use_sensis = false;
        burnin = true;
    case 3
        mName = 'Adaptive MCMC';
        fkt = @mcmc_adaptive;
        adjust_scaling = true;
        use_sensis = false;
        burnin = true;
    case 4
        mName = 'Fisher based MCMC';
        fkt = @mcmc_fish;
        adjust_scaling = false;
        use_sensis = true;
        burnin = false;     
    case 5
        mName = 'Simplified MMALA svdreg';
        fkt = @mcmc_mmala_svd;
        adjust_scaling = false;
        use_sensis = true;
        burnin = false;
    case 6
        mName = 'Simplified MMALA chol';
        fkt = @mcmc_mmala_chol;
        adjust_scaling = false;
        use_sensis = true;
        burnin = false;

end

%% Exchange Method
ExchangeName = 'None';                    
if ( exchange_method == 0 )
    ExchangeName = 'No exchanges';
end
if ( exchange_method == 1 )
    ExchangeName = 'Parallel tempering';
    ex_fkt = @exchange_PT;
end
% if ( exchange_method == 1 )
%     exName = 'Parallel Hierarchical Sampling';
%     cScale = exchange;
%     temp   = ones( 1, n_chains );
%     ex_fkt = @exchange_PHS;
% end



if(~burnin)
    nburnin = 0;
end

jrungo = -nburnin+1;
accept_rate = zeros(NumberOfChains,1);
exchange_rate = zeros(NumberOfChains,1);
% mcmc
arWaitbar(0);
accString                            = repmat( '%4.1f%% ', 1, NumberOfChains );
numString                            = repmat( '%4.1f ', 1, NumberOfChains );

LogPosterior_curr              = nan(NumberOfChains,1);
LogPrior_curr                  = nan(NumberOfChains,1);
LogPosterior_proposal          = nan(NumberOfChains,1);
LogPrior_proposal              = nan(NumberOfChains,1);

TemperedLogPosterior_curr      = nan(NumberOfChains,1);
TemperedLogPosterior_proposal  = nan(NumberOfChains,1);


%accepts                    = nan(NumberOfChains,LengthOfAcceptanceTestChain);
accepts                    = nan(NumberOfChains,nruns*nthinning);
exchanges                  = nan(NumberOfChains,LengthOfAcceptanceTestChain);
i_accepts                  = 1;
para_proposal              = nan(size(para_curr));
mu_curr                     = nan(size(para_curr));
covar_curr                  = nan(length(ar.p(qFitGlobal)),length(ar.p(qFitGlobal)), NumberOfChains);


ProposedSwaps              = zeros(NumberOfChains);
AcceptedSwaps              = zeros(NumberOfChains);


CovarianceAdaptationFactor = zeros(NumberOfChains,1);





if(~use_sensis)
    for chID = 1 : NumberOfChains
        arCalcMerit(false,para_curr(chID,:));
        LogPosterior_curr(chID) = ar.chi2fit;
        LogPrior_curr(chID)     = ar.chi2prior;
        TemperedLogPosterior_curr(chID) = (LogPosterior_curr(chID) - LogPrior_curr(chID))*beta(chID) + LogPrior_curr(chID);                 

        % Calculate tempered distribution value
        %[ LogPosterior_proposal(chID), LogPrior_proposal(chID), ~, ~] = heatResiduals( beta(chID), ar.chi2fit, ar.chi2prior, ResidualType,qFitGlobal );
        [ mu_curr(chID,:), covar_curr(:,:,chID)] = feval(fkt, para_curr(chID,:), Cfactor(chID),accept_rate(chID),max_accept,min_accept,parasHistory(:,:,chID),parasHistory_index, nwindow);

%         arCalcMerit(false,para_curr(chID,:));
%         [ LogPosterior_curr(chID), LogPrior_curr(chID) ] = heatResiduals( beta(chID), ar.chi2fit, ar.chi2prior, ResidualType,qFitGlobal);
%         [ mu_curr(chID,:), covar_curr(:,:,chID) ] = feval(fkt, para_curr(chID,:),Cfactor(chID),accept_rate(chID),max_accept,min_accept,parasHistory(:,:,chID),parasHistory_index, nwindow);
    end
else
    InvProposalPrior     = diag(1./((ar.ub-ar.lb)/2));
    InvProposalPrior     = InvProposalPrior(qFitGlobal,qFitGlobal);
    for chID = 1 : NumberOfChains
        
        arCalcMerit(true,para_curr(chID,:));
        LogPosterior_curr(chID) = ar.chi2fit;
        LogPrior_curr(chID)     = ar.chi2prior;
        TemperedLogPosterior_curr(chID) = (LogPosterior_curr(chID) - LogPrior_curr(chID))*beta(chID) + LogPrior_curr(chID); 
        TemperedRes_curr = ar.res;
        TemperedRes_curr(ResidualType<=2)                    = sqrt(beta(chID)) *ar.res(ResidualType<=2);
        TemperedSRes_curr = ar.sres(:,qFitGlobal);
        TemperedSRes_curr((ResidualType<=2),:)               = sqrt(beta(chID)) *ar.sres((ResidualType<=2),qFitGlobal);

        %[ LogPosterior_proposal(chID), LogPrior_proposal(chID), res_trial, Sres_trial ] = heatResiduals( beta(chID), ar.chi2fit, ar.chi2prior, ResidualType,qFitGlobal, ar.res, ar.sres(:,qFitGlobal) );
        [  mu_curr(chID,:), covar_curr(:,:,chID)] =   feval(fkt, para_curr(chID,:), TemperedRes_curr, TemperedSRes_curr, InvProposalPrior, RegularizationThreshold);       
%         arCalcMerit(true,para_curr(chID,:));
%         [ LogPosterior_curr(chID), LogPrior_curr(chID), res_curr, Sres_curr ] = heatResiduals( beta(chID), ar.chi2fit, ar.chi2prior, ResidualType,qFitGlobal, ar.res, ar.sres(:,qFitGlobal) );
%         [ mu_curr(chID,:), covar_curr(:,:,chID) ] = feval(fkt, para_curr(chID,:), res_curr, Sres_curr, InvProposalPrior, RegularizationThreshold);
    end
end


% Initialise trial variables since blocking will only do partial in-place
% updates
covar_proposal     = covar_curr;
covar_currUSED     = covar_curr;

mu_proposal        = mu_curr;

fprintf( sprintf( 'Sampling Method: %s\n', mName ) );
fprintf( sprintf( 'Exchange Method: %s\n', ExchangeName ) );
% fprintf( sprintf( 'Scales: %s\n', numString ), cScale );
fprintf( sprintf( 'Inverse Temperatures: %s\n', numString ), beta );
fprintf('MCMC sampling...\n');
tic;
count_chain_reset = 0;
jthin = 1;
jcount = 1;
for jruns = 1:floor(((nruns*nthinning)+nburnin))
    
    % For each chain
    for chID = 1 : NumberOfChains
        para_proposal(chID,:)          = para_curr(chID,:);
        mu_proposal(chID,:)            = mu_curr(chID,:);
        covar_proposal(:,:,chID)       = covar_curr(:,:,chID);
        
        % Calculate and update Acceptance Rate
%         accept_rate(chID) = sum(accepts(chID,:))/LengthOfAcceptanceTestChain;
        accept_rate(chID) = nansum(accepts(chID,:))/sum(~isnan(accepts(chID,:)));
        if ( exchange_method > 0 )
            exchange_rate(chID) = sum(exchanges(chID,:))/LengthOfAcceptanceTestChain;
        end
        if(chID==1)
            if(jrungo>0)
                arWaitbar(jruns, floor((nruns*nthinning)+nburnin), sprintf( sprintf('MCMC run (acceptance rates %s)', accString), ...
                    accept_rate*100));
            else
                arWaitbar(jruns, floor((nruns*nthinning)+nburnin), sprintf( sprintf('MCMC burn-in (acceptance rates %s)', accString), ...
                    accept_rate*100));
            end
        end
        
        % Calculate Proposal    
        para_proposal(chID,locs)        = mvnrnd(mu_curr(chID,locs), ManualScalingFactor*covar_currUSED(locs,locs,chID));
        LogPosterior_proposal(chID)     = inf;

        
        % Check Bounds
        if((sum(para_proposal(chID,locs)<lb(locs)) + sum(para_proposal(chID,locs)>ub(locs))) == 0) 
            % Calculate LogLikelihood and LogPrior
            try
                if(~use_sensis)
                    arCalcMerit(false,para_proposal(chID,:));
                    LogPosterior_proposal(chID) = ar.chi2fit;
                    LogPrior_proposal(chID)     = ar.chi2prior;
                    TemperedLogPosterior_proposal(chID) = (LogPosterior_proposal(chID) - LogPrior_proposal(chID))*beta(chID) + LogPrior_proposal(chID);                 

                    % Calculate tempered distribution value
                    %[ LogPosterior_proposal(chID), LogPrior_proposal(chID), ~, ~] = heatResiduals( beta(chID), ar.chi2fit, ar.chi2prior, ResidualType,qFitGlobal );
                    [ mu_proposal(chID,locs), covar_proposal(locs,locs,chID)] = feval(fkt, para_proposal(chID,locs), Cfactor(chID),accept_rate(chID),max_accept,min_accept,parasHistory(:,:,chID),parasHistory_index, nwindow);

                else
                    arCalcMerit(true,para_proposal(chID,:));
                    LogPosterior_proposal(chID) = ar.chi2fit;
                    LogPrior_proposal(chID)     = ar.chi2prior;
                    TemperedLogPosterior_proposal(chID) = (LogPosterior_proposal(chID) - LogPrior_proposal(chID))*beta(chID) + LogPrior_proposal(chID); 
                    TemperedRes_proposal = ar.res;
                    TemperedRes_proposal(ResidualType<=2)                    = sqrt(beta(chID)) *ar.res(ResidualType<=2);
                    TemperedSRes_proposal = ar.sres(:,qFitGlobal);
                    TemperedSRes_proposal((ResidualType<=2),:)               = sqrt(beta(chID)) *ar.sres((ResidualType<=2),qFitGlobal);

                    %[ LogPosterior_proposal(chID), LogPrior_proposal(chID), res_trial, Sres_trial ] = heatResiduals( beta(chID), ar.chi2fit, ar.chi2prior, ResidualType,qFitGlobal, ar.res, ar.sres(:,qFitGlobal) );
          
                    
                    [ mu_proposal(chID,locs), covar_proposal(locs,locs,chID)] =   feval(fkt, para_proposal(chID,locs), TemperedRes_proposal, TemperedSRes_proposal, InvProposalPrior, RegularizationThreshold);
                end
                
                ProbProposalGivenCurr(chID)  = mvnpdf(para_proposal(chID,locs), mu_curr(chID,locs), covar_curr(locs,locs,chID));
                ProbCurrGivenProposal(chID)  = mvnpdf(para_curr(chID,locs), mu_proposal(chID,locs), covar_proposal(locs,locs,chID));
                a = exp(-0.5*(TemperedLogPosterior_proposal(chID) - TemperedLogPosterior_curr(chID))) * (ProbCurrGivenProposal(chID) / ProbProposalGivenCurr(chID));
                
                randa = rand;
                qa = randa <= min([1 a]); % accept?

            catch
                LogPosterior_proposal(chID)           = inf;
                LogPrior_proposal(chID)               = inf;
                TemperedLogPosterior_curr(chID)       = inf;
                TemperedLogPosterior_proposal(chID)   = inf;
                
                disp( 'Simulation Failed => Rejecting step' );
                qa = 0;
            end
        else
            %fprintf('#%i bound violation\n', jrungo);
            qa = false;
            NumberOfBoundViolations = NumberOfBoundViolations+1;
            %BoundViolationIterations(NumberOfBoundViolations)=jrungo;
            %L_trial(chID) = inf;
            %P_trial(chID) = inf;
        end


        
        
        %% update structs and variables
        if(qa)
            para_curr(chID,:)                = para_proposal(chID,:);
            LogPosterior_curr(chID)          = LogPosterior_proposal(chID);
            TemperedLogPosterior_curr(chID)  = TemperedLogPosterior_proposal(chID);
            LogPrior_curr(chID)              = LogPrior_proposal(chID);
            mu_curr(chID,:)                  = mu_proposal(chID,:);
            covar_curr(:,:,chID)             = covar_proposal(:,:,chID);

            if(method==3)
                parasHistory(parasHistory_index,:,chID) = para_curr(chID,:);
                if ( chID == NumberOfChains )
                    parasHistory_index = parasHistory_index + 1;
                    if(parasHistory_index > nwindow)
                        parasHistory_index = 1;
                    end
                end
            end
        else

        end
        
        accepts(chID,i_accepts) = qa;
        

                %% Scaling Adjustment        
        if UseScaling == 1
            % Scaling adjustment from Lacki, Miasojedow: State-dependent swap
            % strategies and atomatic reduction of number of temperatures in
            % adaptive parallel tempering algorithm

            % Using stationary covariance (Haario et al. 2001)
            % Initialize parameters:
            %MemoryLength = 1;  % Number of steps for which almost no adaptation is made in the beginning
             % Control parameter for decay of Covariance Adaptation, (0<p<1), higher values mean fast decay 
%             i = max(jruns+1,MemoryLength);
% 
%             Gamma = 1/(i^DecayParameter);
% 
%             % Adaptation
%             mu_curr(chID,:) = (1-Gamma)*mu_curr(chID,:) + Gamma*para_curr(chID,:);
%             covar_curr(:,:,chID)  = (1-Gamma)*covar_curr(:,:,chID)+ Gamma*(para_curr(chID,:)-mu_curr(chID,:))*(para_curr(chID,:)-mu_curr(chID,:))';         
           if isnan(accept_rate(chID))
           else
             CovarianceAdaptationFactor(chID) = CovarianceAdaptationFactor(chID)+((accept_rate(chID))-0.234)*((jruns+1)^(-DecayParameter));
             % Adjust covariance matrix (propsed and current
             covar_currUSED(:,:,chID) = exp(CovarianceAdaptationFactor(chID))*covar_curr(:,:,chID);

             % Regularization of proposal covariance
%              feval(fkt, para_curr(chID,:), Cfactor(chID),accept_rate(chID),max_accept,min_accept,parasHistory(:,:,chID),parasHistory_index, nwindow);
%              mcmc_mmala_chol(ptmp, restmp, srestmp,InvProposalPriorTemp, RegularizationThresholdTemp)

             %ParameterNumber = length(para_curr(chID,:));

%             [~,qq] = cholcov(covar_currUSED(:,:,chID),0);
%             if qq ~= 0
%                covar_currUSED(:,:,chID) = covar_currUSED(:,:,chID) + RegularizationThreshold*eye(ParameterNumber);
%                covar_currUSED(:,:,chID) = (covar_currUSED(:,:,chID)+covar_currUSED(:,:,chID)')/2;
%                [~,qq] = cholcov(covar_currUSED(:,:,chID),0);
%                if qq ~= 0
%                   covar_currUSED(:,:,chID) = covar_currUSED(:,:,chID) + max(max(covar_currUSED(:,:,chID)))/1000*eye(ParameterNumber);
%                   covar_currUSED(:,:,chID) = (covar_currUSED(:,:,chID)+covar_currUSED(:,:,chID)')/2;
%                end
%             end
           end
        else
            covar_currUSED(:,:,chID) = covar_curr(:,:,chID);
        end                

        
        
        
        % adjust scaling during burn-in (OLD SCALING ADJUSTMENT)
        if(jrungo<=0 && adjust_scaling)
            if(accept_rate(chID) > max_accept && Cfactor(chID)*Cmod<Cmax)
                Cfactor(chID) = Cfactor(chID)*Cmod;
            elseif(accept_rate(chID) < min_accept && Cfactor(chID)/Cmod>Cmin)
                Cfactor(chID) = Cfactor(chID)/Cmod;
            end
        end
        
        
        

    end
    
    i_accepts = i_accepts + 1;
    if(i_accepts>length(accepts))
        i_accepts = 1;
    end 

%     % Save Samples
%     if(jrungo>0)
%         jthin = jthin + 1;
%         if(jthin > nthinning)
%             ar.ps(jcount+jindexoffset,:,:) = paraReset;
%             ar.ps(jcount+jindexoffset,qFitGlobal,:) = para_curr.';
%             ar.ps_trial(jcount+jindexoffset,:,:) = paraReset;
%             ar.ps_trial(jcount+jindexoffset,qFitGlobal,:) = para_proposal.';
%             ar.chi2s(jcount+jindexoffset,:) = LogPosterior_curr;
%             ar.chi2s_trial(jcount+jindexoffset,:) = LogPosterior_proposal;
%             ar.Q_trial(jcount+jindexoffset,:) = ProbProposalGivenCurr;
%             ar.Q_curr(jcount+jindexoffset,:)  = ProbCurrGivenProposal;
%             ar.acceptance(jcount+jindexoffset,:) = accept_rate;
%             if ( exchange_method > 0 )
%                 ar.exchange(jcount+jindexoffset,:) = exchange_rate;
%             end
%             jthin = 1;
%             jcount = jcount + 1;
%         end
%     end
%     jrungo = jrungo + 1;
    
    
    
    
    % Perform chain exchange
    if ( exchange_method > 0 )
        if ( NumberOfChains > 1 )
            
            
            SwapProbForward = zeros(length(LogPosterior_curr));
            
            for index1 = 2:NumberOfChains
                for index2 = 1:index1-1
                    SwapProbForward(index2,index1) = exp(-abs(LogPosterior_curr(index2)-LogPosterior_curr(index1)));
                end
            end
            if sum(SwapProbForward(:)) ~= 0
                SwapProbForward = SwapProbForward/sum(SwapProbForward(:));   
            else
                SwapProbForward = ones(length(LogPosterior_curr));
                SwapProbForward = SwapProbForward/sum(SwapProbForward(:));   
            end
            
            SwapIndex = find(cumsum(SwapProbForward(:)) > rand(), 1, 'first');

            %SwapProbBackward = SwapProbForward;

            [k2,k1] = meshgrid(1:NumberOfChains,1:NumberOfChains);   
            k1 = k1(SwapIndex);
            k2 = k2(SwapIndex);
            ProbAccSwap = exp(-0.5*(LogPosterior_curr(k1)-LogPosterior_curr(k2)))^(beta(k2)-beta(k1));
            
            
            % Update chain states and run statistics
            ProposedSwaps(k1,k2) = ProposedSwaps(k1,k2) + 1;
            if rand <= ProbAccSwap
               AcceptedSwaps(k1,k2)   = AcceptedSwaps(k1,k2) + 1;                          
               para_curr([k1,k2],:) = para_curr([k2,k1],:);
               LogPosterior_curr([k1,k2]) = LogPosterior_curr([k2,k1]);
            end
      

%             [ exaccept ] = feval( ex_fkt, fkt );
%             exchanges(:,i_accepts) = exaccept;
%             save_samples( );
%             jrungo = jrungo + 1;
        end
    end
    
    % Save Samples
    if(jrungo>0)
        jthin = jthin + 1;
        if(jthin > nthinning)
            ar.ps(jcount+jindexoffset,:,:) = paraReset;
            ar.ps(jcount+jindexoffset,qFitGlobal,:) = para_curr.';
            ar.ps_trial(jcount+jindexoffset,:,:) = paraReset;
            ar.ps_trial(jcount+jindexoffset,qFitGlobal,:) = para_proposal.';
            ar.chi2s(jcount+jindexoffset,:) = LogPosterior_curr;
            ar.chi2s_trial(jcount+jindexoffset,:) = LogPosterior_proposal;
            ar.Q_trial(jcount+jindexoffset,:) = ProbProposalGivenCurr;
            ar.Q_curr(jcount+jindexoffset,:)  = ProbCurrGivenProposal;
            ar.acceptance(jcount+jindexoffset,:) = accept_rate;
            if ( exchange_method > 0 )
                ar.exchange(jcount+jindexoffset,:) = exchange_rate;
                ar.AcceptedSwaps = AcceptedSwaps;
                ar.ProposedSwaps = ProposedSwaps;
            end
            jthin = 1;
            jcount = jcount + 1;
        end
    end
    jrungo = jrungo + 1;   
    
    
end

if ( fitErr > 0 )
    disp( 'Set UseFitCorrection back to its original value' );
    ar.useFitErrorCorrection = fitErr;
end

arWaitbar(-1);
fprintf('done (%s, %i chain resets) \n', secToHMS(toc), count_chain_reset);
ar.mc3.toc = toc;


mcmc.ps         = ar.ps;
mcmc.ps_trial   = ar.ps_trial;
mcmc.chi2s      = ar.chi2s;
mcmc.chi2s_trial= ar.chi2s_trial;
mcmc.acceptance = ar.acceptance;
mcmc.MethodName = mName;
if ( exchange_method > 0 )
mcmc.exchange   = ar.exchange;
mcmc.AcceptedSwaps = ar.AcceptedSwaps;
mcmc.ProposedSwaps = ar.ProposedSwaps;
end

fprintf('\n%8.0f bound violations observed (for all chains)\nThis is %3.0f percent of all iterations.\n\n', NumberOfBoundViolations,100*NumberOfBoundViolations/(nruns*NumberOfChains));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sampling functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% N(0,1)
    function [mu, covar] = mcmc_norm_mvnrnd(ptmp ,~,~,~,~,~,~,~)
        mu = ptmp;
        covar = eye(length(ptmp));
    end
    
% N(0,c) scaled
    function [mu, covar] = mcmc_scaled_mvnrnd(ptmp, CfactorTemporary,~,~,~,~,~,~)
        mu = ptmp;
        covar = eye(length(ptmp)) * CfactorTemporary;
    end
    
% Adaptive MCMC
    function [mu, covar] = mcmc_adaptive(ptmp, CfactorTempo,AcceptRateTemp,max_accept,min_accept,ps_histTemp,ps_hist_index, nwindow)
%         fprintf('%i/%i %e\n', sum(isnan(ps_hist(:,1))), nwindow, Cfactor);
       if((isnan(AcceptRateTemp) || AcceptRateTemp > max_accept || AcceptRateTemp < min_accept) ...
                && (sum(isnan(ps_histTemp(:,1)))==0))% || jrungo<=0))
            mu = ptmp;            
            if(sum(isnan(ps_histTemp(:,1))) == 0)
                % muAdapt = mean(ps_hist);
                covar = cov(ps_histTemp(:,:))* CfactorTempo;
            else
                covar = eye(length(ptmp)) * CfactorTempo;
            end            
        else
            if(ps_hist_index == nwindow && sum(isnan(ps_histTemp(:,1)))==0)
                % muAdapt = mean(ps_hist);
                covar = cov(ps_histTemp(:,:))* CfactorTempo;
%                 figure(1)
%                 plot(ps_hist(:,1),'x-');
%                 figure(2)
%                 imagesc(CAdapt);
%                 colorbar
            else
                covar = eye(length(ptmp)) * CfactorTempo;
            end
            mu = ptmp;
            % mu = muAdapt;
%             covar = CAdapt(:,:,chID);
        end
    end
    
    
        
% (Simple) Fisher based
    function [mu, covar] = mcmc_fish(ptmp, ~, srestmp, InvProposalPriorTemp, ~)
        GMetric = srestmp'*srestmp;
        alpha_dash = GMetric + InvProposalPriorTemp;
        
        mu=ptmp;
        covar=inv(alpha_dash);
    end

% MMALA (simplified)Cholesky regularization
    function [mu, covar] = mcmc_mmala_chol(ptmp, restmp, srestmp,InvProposalPriorTemp, RegularizationThresholdTemp)
        Gradient = - restmp*srestmp;
        GMetric = srestmp'*srestmp;
        alpha_dash = GMetric + InvProposalPriorTemp;
        
        ParameterNumber = length(ptmp);
        
       % Regularize the Metric Tensor
       [~,p] = cholcov(alpha_dash,0);
       k = 0;
       if p ~= 0          
          while p ~= 0
             alpha_dash_k = alpha_dash + 10^k*RegularizationThresholdTemp*eye(ParameterNumber);
             [~,p] = cholcov(alpha_dash_k,0);
             k = k+1;
          end
          alpha_dash = alpha_dash_k;     
       end
        
        mu = ptmp + (1/2)*mldivide(alpha_dash,Gradient')';
        covar = inv(alpha_dash);
        covar = 0.5*(covar+covar');
        
        % Regularization of proposal covariance
        [~,p] = cholcov(covar,0);
        if p ~= 0
           covar = covar + RegularizationThresholdTemp*eye(ParameterNumber);
           covar = (covar+covar')/2;
           [~,p] = cholcov(covar,0);
           if p ~= 0
              covar = covar + max(max(covar))/1000*eye(ParameterNumber);
              covar = (covar+covar')/2;
           end
        end
        
end
    
    
    
% MMALA (simplified) SVD Regularization
    function [mu, covar] = mcmc_mmala_svd(ptmp, restmp, srestmp,InvProposalPriorTemp, RegularizationThresholdTemp)
        Gradient    = - restmp*srestmp;
        GMetric     = srestmp'*srestmp;
        alpha_dash  = GMetric + InvProposalPriorTemp;

        % solve with SVD regularization
        [U,S,V] = svd(alpha_dash);
        s = diag(S);
        qs = (s/max(s)) > RegularizationThresholdTemp;

        deltap = zeros(size(Gradient));
        for jj=find(qs)'
            deltap = deltap + transpose((U(:,jj)'*Gradient'/S(jj,jj))*V(:,jj));
        end
        mu      = ptmp + 0.1*(1/2)*deltap;
        covar   = inv(alpha_dash);              
    end
    

