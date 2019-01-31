%ARREGULARIZE The major analysis function for regularization of a given system.
% 
% To obtain reltos of dynamic parameters between e.g. two cell lines, run
% l1_PrintDynPars. The output can be parsed into the CONDITIONS sectino of
% the model.def file
%
% [runtime,parsimonious_model] = ARREGULARIZE(RELTOS,LAMBDAMIN,LAMBAMAX, VARARGIN)
%
%   reltos      -either char that is included exactly in all fold-change
%               labels (ar.pLabel), e.g. 'relto' or 'fc'
%               -or a vector indicating the indices of ar.p that should be
%               penalized
%   lambdamin   float, lower bound of penalty strength scan range. NB: in
%               logarithmic lambda space, this is the order of magnitude
%               rather than the actual number, e.g. '1' for lambda_min=10
%               in log-space or lambda_min=1 in linear space, resp.
%   lambdamax   float, upper bound of penalty strength scan range. NB: see
%               lambdamin
%
%   lambdanum   integer: positive number indicates the total number of
%   (optional)  lambda strengths to be scanned. Negative numbers indicate
%               the number of steps per order of magnitude (log scale) or
%               per one-step in lambda. Default: -6
%               Ex.: arRegularize('relto',0,4) penalizes from 1 to 10000 in
%               4 * 6 + 1 = 25 log steps, arRegularize('relto',0,4,25) does
%               the same.
% 
% varargin: Name,Value options, e.g.:
%       'LambdaScale': {'lin','log'}, determines whether the penalty strength
%                   are chosen as linspace(lambdamin, lambdamax) or as
%                   logspace(lambdamin,lambdamax). Default: 'log'
% 
%       'IncludePriored': true or false (def.): determines whether other
%                   parameter priors for fold-changes (!) are overwridden
% 
%       'ReleaseFixed': true or false(def.): determines whether parameters that
%                   are fixed for fits should be released to be fitted during
%                   regularization
% 
%       'Type'      'l1' (Def.) for original l1-penalization,
%                   'lq' for lq penalization
%                   'al', 'adaptive', 'adaptive-GW' for Adaptive (Group) Lasso
%                       penalty: 'adaptive-GW' uses group-wise weights, whereas
%                       'adaptive' and 'al' yield parameter-wise weights. For
%                       ungrouped regularization, there is no difference
%                   'en' for Elastic Net
% 
%       'Deform'    Float, Def. 0.2: Deviation of standard l1, indicates
%                   adaptivity gamma, elasticity alpha or 1-q for lq penalties.
%                   'Type','weighted' sets Deform to -1, so the weights
%                   provided as 'Weights' are conserved
% 
%       'Sequential' {'none','sorted',lexicographic'}. Default: 'none'. For
%                   simultaneous scans, provide 'none', for sequential scan
%                   with parameters tested against zero with ascending modulus,
%                   try 'sorted', for sequential scan with order as is, use
%                   'lexicographic'.
% 
%       'Weights'   'OLS' or numeric array. Default: 1. If 'Type' is provided 
%                   as Adaptive(Group) Lasso, then 'Weights' is set to 'OLS'
% 
%       'GroupAlign' numeric vector or string. Default: 'none'
%                   vector: should have the length of the relto-vector,
%                       indicates group numbers for each index
%                       Ex.: [1 2 1 1 2] -> parameters with indices            
%                       l1pars(1), l1pars(3), l1pars(4) belong to the first
%                       group. parameters with indices l1pars(2) and l1pars(5)
%                       belong to the second group.
%                   'altX', 'groupX','single' assign group numbers as follows:
%                       for 'altX', every X-th parameter belongs to the same
%                       group; for 'groupX', X subsequent parameters are
%                       grouped together; 'single' provides one-parameter
%                       groups (actually original l1, it's better not to use
%                       this option except for testing purposes)
%                       Ex.: 'alt3' -> groupalign = [1 2 3 1 2 3 1 2 3 ...],
%                       'group3' -> groupalign = [1 1 1 2 2 2 3 3 3 ...],
%                       'single' -> groupalign = 1:length(l1pars)
%                   'n:m' with n,m integer or end-integer. Then, the group 
%                       labels (ar.pLabel) are split at the character '_'. 
%                       Then the substrings with indices n:m are collected an
%                       parameters having the same results are grouped
%                       together.
%                       Example: ar.pLabel =
%                       {'relto_A_0_CL_1', 'relto_A_0_CL_2', 
%                        'relto_k_1_CL_1', 'relto_k_1_CL_2',
%                        'relto_k_2_CL_1', 'relto_k_2_CL_2'};
%                        'GroupAlign','2:3' OR 'GroupAlign','alt2' groups by
%                        the strings 'A_0','k_1','k_2' (3 groups),
%                        whereas 'GroupAlign','4:5' OR
%                        'GroupAlign','end-1:end', OR 'GroupAlign','group2'
%                        would group by 'CL_1','CL_2' (2 groups).
% 
%       'Test'      {'LRT','BIC'}, Def. 'LRT'. Determines whether Likelihood
%                   Ratio Test or Bayesian Information Criterion are used to
%                   find the parsimonious model
% 
%       'Means'     array compatible with relto-vector length, or 'keep'.
%                   Default: 0. Sets the target value/mean value of the
%                   penalization. 'keep' does not change a mean value provided
%                   beforehand
% 
%       'UpperBounds', 'LowerBounds': 'keep' (def.) or array compatible with
%                   relto-vector length. Sets parameter constraints for
%                   fold-changes. Behaves like 'Means'
% 
%       'Gradient'  float. Default: 0. Adds a gradient to the optimization, see
%                   l1 functions for further explanation.
% 
%       'Threshold' positive float. Default: 1e-6. Determines epsilon > 0 below
%                   which a fold-change is considered as zero when shrinking
%                   the model
% 
%       'OptimizerSteps' array containing at position i, how many steps the
%                   optimization with ar.config.optimizer = i shall do at most.
%                   Default [1000 20] (1000 steps lsqnonlin, 20 steps fmincon)
% 
%       'doPlot','doTree': true, false. Def.: 'doPlot',true,'doTree',false.
%                   Determines whether l1/grplasPlot and/or l1/grplasTree are
%                   called
% 
%   runtime     nan: if the regularization fails
%   (optional)  float > 0 if regularization successful, indicating the run 
%               time in seconds
%
%   parsimonious_model: indices of the final model which are left free
%   (optional)  for optimization (ar.qFit == 1)
%
% 
% See also l1Init, l1Scan, l1Seq

function  varargout  = arRegularize( reltos, lambdamin, lambdamax, varargin )
switch nargout
    case 0
        varargout = {};
    case 1
        varargout = {nan};
    case 2
        varargout = {nan,[]};
end

global ar;

if(isempty(ar))
    fprintf('please initialize by arInit')
    return
end

%% Parsing input
p = inputParser;

validType = {'l1','lq','adaptive','elasticnet','AL','EN',...
    'adaptive-GW','weighted'};
validLambdaScale = {'log','lin'};
validSequential = {'none','sorted','lexicographic'};
validTest = {'LRT','BIC'};

checkNumericVector = @(x) isnumeric(x) && isvector(x);
checkNumericScalar = @(x) isnumeric(x) && isscalar(x);

checkType = @(x) any(validatestring(x,validType)) || ismatrix(x);
checkReltos = @(x) ischar(x) || checkNumericVector(x);
checkKeep = @(x) strcmpi(x,'keep') || checkNumericVector(x);
checkLarger0 = @(x) checkNumericScalar(x) && x>0;
checkDeform = @(x) all(isnumeric(x));
checkWeights = @(x) strcmpi(x,'OLS') || all(isnumeric(x) & x>0);
checkGroupAlign = @(x) ischar(x) || ...
    checkNumericVector(x) || iscell(x);
checkLambdaScale = @(x) any(validatestring(x,validLambdaScale));
checkSequential = @(x) any(validatestring(x,validSequential));
checkTest = @(x) any(validatestring(x,validTest));

addRequired(p,'Reltos',checkReltos)
addRequired(p,'LambdaMin',checkNumericScalar)
addRequired(p,'LambdaMax',checkNumericScalar)
addOptional(p,'LambdaNum',-6,checkLarger0)

addParameter(p,'LambdaScale','log',checkLambdaScale)
addParameter(p,'Type','l1',checkType)
addParameter(p,'IncludePriored',false,@islogical)
addParameter(p,'ReleaseFixed',false,@islogical)
addParameter(p,'Means',0,checkKeep)
addParameter(p,'UpperBounds','keep',checkKeep)
addParameter(p,'LowerBounds','keep',checkKeep)
addParameter(p,'Threshold',1e-6,checkLarger0)
addParameter(p,'Deform',nan,checkDeform)
addParameter(p,'Weights',1,checkWeights)
addParameter(p,'GroupAlign','none',checkGroupAlign)
addParameter(p,'OptimizerSteps',[1000 20],checkNumericVector)
addParameter(p,'Sequential','none',checkSequential)
addParameter(p,'Test','LRT',checkTest)
addParameter(p,'Gradient',0,checkNumericScalar)
addParameter(p,'doPlot',true,@islogical)
addParameter(p,'doTree',false,@islogical)

parse(p, reltos, lambdamin, lambdamax, varargin{:})

type = validatestring(p.Results.Type,validType);
groupalign = p.Results.GroupAlign;
test = validatestring(p.Results.Test,validTest);
lambdascale = validatestring(p.Results.LambdaScale,validLambdaScale);
sequential = validatestring(p.Results.Sequential,validSequential);
weights_user = p.Results.Weights;
deform = p.Results.Deform;

if (strcmpi(type,'al') || strcmpi(type,'adaptive') || strcmpi(type,...
        'adaptive-GW'))
    weights_user = 'OLS';
elseif strcmpi(type,'weighted')
    deform = -1;
    type = 'AL';
end

%% Determining the range of penalty strength
try
    
    lammin = p.Results.LambdaMin;
    lammax = p.Results.LambdaMax;
    if lammax < lammin
        error('LambdaMax < LambdaMin, so scan range is empty')
    end
    if p.Results.LambdaNum > 0
        nticks = p.Results.LambdaNum;
        % User-supplied number of penalty strengths
    else
        if strcmpi(p.Results.LambdaScale,'lin')
            nticks = ceil((log10(lammax) - log10(lammin)) * (-p.Results.LambdaNum)) + 1;
        else
            nticks = ceil((lammax - lammin) * (-p.Results.LambdaNum)) + 1;
        end
        % User-supplied number of penalty strengths per order of magnitude
    end
    if nticks > 1
        if strcmpi(lambdascale,'log')
            linv = logspace(lammin,lammax,nticks);
            % NB: lammin, lammax refer to the order of magnitude
        elseif strcmpi(lambdascale,'lin')
            linv = linspace(lammin,lammax,nticks);
        else
            error('Unrecognized Lambda Scale')
        end
    elseif (nticks == 1 && lammax == lammin)
        linv = lammax;
        % Case of only one lambda value (less useful anyway)
    else
        error('Unclear Lambda range detected: lambdanum = 1, but lambdamax <> lambdamin')
    end

    linv = 1./linv;
    % Penalty strength must be passed to scan routines as 1/lambda
catch ME
    fprintf('Lambda Values could not be assigned:\n')
    fprintf('%s\n',ME.message)
    return
end

%% Determining parameters to be penalized
try
    if ischar(p.Results.Reltos)
        l1pars = arPrint(p.Results.Reltos);
        % Returns all parameters with labels containing p.Results.Reltos
    elseif isvector(p.Results.Reltos)
        l1pars = 1:length(ar.p);
        l1pars = l1pars(p.Results.Reltos);
    else
        error('Malformed fold-change specification, please provide chars or num. array')
    end
    fixed = ar.qFit(l1pars) == 2;
    if any(fixed)
        if p.Results.ReleaseFixed
            ar.qFit(l1pars) = 1;
            fprintf('Reset parameters %s, refitting the model! \n',mat2str(l1pars(fixed)))
            arFit
        else
            warning('The following indices are not fitted (ar.qFit = 2): %s.\n If you would like to penalize them anyway, use "ReleaseFixed", true\n',...
                mat2str(l1pars(fixed)))
            l1pars = l1pars(~fixed);
        end
    end
    priored = ar.type(l1pars) ~= 0;
    if any(priored)
        if p.Results.IncludePriored
            ar.type(l1pars) = 0;
            fprintf('Reset priors for parameters %s, refitting the model! \n',mat2str(l1pars(~priored)))
            arFit
        else
            warning('The following indices are already subject to priors (ar.type <> 0): %s.\nIf you would like to penalize them anyway, use "IncludePriors", true\n',...
                mat2str(l1pars(priored)))
            l1pars = l1pars(~priored);
        end
    end
catch ME
    fprintf('Please use a valid fold-change identifier:\n')
    fprintf('%s\n',ME.message)
    return
end

%% Determining Means, Upper and Lower Bonds
try
    % 'keep' (def.) takes the existing properties of ar with respect to
    % means, upper and lower bounds
    
    if strcmpi(p.Results.Means,'keep')
        means = ar.means(l1pars);
    else
        means = p.Results.Means .* ones(size(l1pars));
    end
    if strcmpi(p.Results.LowerBounds,'keep')
        lbs = ar.lb(l1pars);
    else
        lbs = p.Results.LowerBounds .* ones(size(l1pars));
    end
    if strcmpi(p.Results.UpperBounds,'keep')
        ubs = ar.ub(l1pars);
    else
        ubs = p.Results.UpperBounds .* ones(size(l1pars));
    end
catch ME
    fprintf('Means, upper or lower bounds not compatible with fold-changes:\n')
    fprintf('%s\n',ME.message)
    return
end

if strcmpi(groupalign,'none')
    % Parameter-Wise
    
    %% Determining deformation parameter and penalty function type
    
    try
        weights = 1;
        if strcmpi(type,'lq')
            if isnan(deform)
                metapar = 0.8;
                % heuristically found to be a reasonable choice
            else
                metapar = 1 - deform;
                % Deform is 1 - Exponent here!
            end
            if (metapar <= 0  || metapar >= 1)
                error('Incompatible Exponent, use 0 < q < 1')
            end
        elseif (strcmpi(type,'AL') ...
                || strcmpi(type,'adaptive')...
                || strcmpi(type,'adaptive-GW'))
            if isnan(deform)
                metapar = 0.2;
                % heuristically found to be a reasonable choice
            else
                metapar = deform;
            end
            if strcmpi(weights_user, 'OLS')
                try
                    arFit(true)
                    weights = ar.p(l1pars);
                catch ME
                    fprintf('OLS Weight assignment not possible, fit failed:\n')
                    fprintf('%s\n',ME.message)
                    return
                end
            else
                weights = weights_user;
                % User supplied weights
            end
        elseif (strcmpi(type,'EN') || strcmpi(type,'elasticnet'))
            if isnan(deform)
                metapar = 0.05;
                % possibly a good choice (depends on model)
            else
                metapar = deform;
            end
            if (metapar > 1 || metapar < 0 )
                error('Incompatible Elasticity, use 0 <= alpha <= 1')
            end
        else
            metapar = 0;
        end
    catch ME
        fprintf('The deformation parameter was rejected:\n')
        fprintf('%s\n',ME.message)
        return
    end
    
    %% Actual Scan
    try
        fprintf('Starting %s-Regularization,\n',type)
        fprintf('Penalizing %i fold-changes,\n',length(l1pars))
        if strcmpi(lambdascale,'log')
            addstr = '10^';
        else
            addstr = '';
        end
        fprintf('Scanning over %i penalty strengths from %s%.2g to %s%.2g\n',...
            nticks,addstr,lammin,addstr,lammax)
        
        tic
        l1Init(l1pars,means,lbs,ubs,linv,p.Results.Threshold,type,metapar,weights,false)
        
        if strcmpi(sequential,'none')
            fprintf('in simultaneous penalization mode\n')
            l1Scan(l1pars,[],p.Results.Gradient,[],p.Results.OptimizerSteps)
        elseif strcmpi(sequential,'sorted')
            fprintf('in sequential penalization mode, ordered by absolute value\n')
            l1Seq(l1pars,[],p.Results.Gradient,[],1)
        elseif strcmpi(sequential,'lexicographic')
            fprintf('in sequential penalization mode, ordered as is\n')
            l1Seq(l1pars,[],p.Results.Gradient,[],2)
        end
        
        l1Unpen(l1pars)
        l1SelectOpt(l1pars,test)
        tmp_scantime = toc;
        
        fprintf('Regularization was successful in %.3f sec.\n',tmp_scantime)
        if ar.L1final_ind == 1
            fprintf('First penalty strength led to parsimonious model\n')
            fprintf('Consider decreasing lambdamin\n')
        end
        if ar.L1final_ind == length(linv)
            fprintf('Last penalty strength led to parsimonious model\n')
            fprintf('Consider increasing lambdamax\n')
        end
        
        if p.Results.doPlot
            l1Plot
        end
        if p.Results.doTree
            l1Tree
        end
    catch ME
		fprintf('Scan routine for parameter-wise scan failed:\n')
        fprintf('%s\n',ME.message)
        return
    end
    
else
    % Group Lasso Scans
    
    %% Determine Grouping
    
    % Group-Alignment specification is a bit tricky. The parameter labels
    % are (in mente) subdivided into blocks, which are separated by _. As a
    % string/char array, the blocks that should be considered for grouping
    % can be provided by '1', '2:end-1', 'end' etc. Bracket notation such 
    % as '[1 2 3]' is not supported yet. This nomenclature has 
    % to be consistent among all labels. If the scheme is
    % 'relto_parametername_cellline', then providing '2' groups by
    % parameter if parameternames do not contain underscores. '2:end-1'
    % would do the same job if cell line specifications do not contain _. 
    % If both contain '_', consider renaming them, or proceed to this
    % alternative: It is
    % possible to pass a vector of length(l1pars) which encodes the group
    % assignment for each fold-change to penalize. In addition, it is 
	% possible to use 'altX' or 'groupX' to group every X-th parameter or
	% the X subsequent parameters, respectively.
    
    try
        l1labels = ar.pLabel(l1pars);
        if ischar(groupalign) 
            if (~contains(groupalign,'alt') && ...
                    ~contains(groupalign,'group') && ...
                    ~contains(groupalign,'single'))
                margs = strsplit(groupalign,':');
                % Split up if a range is provided.
                if length(margs) == 1
                    % No range, only a single value
                    ind = str2double(margs{1});
                    if isnan(ind)
                        error('Group Align input index not recognized: Conversion to numeric value of "%s" failed',...
                            margs{1})
                    else
                        inds = [ind ind];
                    end
                    % double storage necessary for later calculations
                elseif length(margs) == 2
                    % A possible range indication
                    
                    tmpl = strsplit(margs{1},'-');
                    % Starting point of range, extract 'end'
                    
                    tmpr = strsplit(margs{2},'-');
                    % End point of range, extract 'end'
                    
                    if length(tmpl) == 1
                        if contains(tmpl{1},'end')
                            % 'end' was provided
                            indl = 0;
                        else
                            % No 'end', no subtraction
                            indl = str2double(tmpl{1});
                        end
                    elseif (length(tmpl) == 2 && contains(tmpl{1},'end'))
                        % Subtraction with first term 'end'
                        indl = str2double(tmpl{2});
                        indl = -indl;
                    else
                        error('Group Align: Left Boundary malformed')
                    end
                    
                    if length(tmpr) == 1
                        if contains(tmpr{1},'end')
                            % 'end' was provided
                            indr = 0;
                        else
                            % No 'end', no subtraction
                            indr = str2double(tmpr{1});
                        end
                    elseif (length(tmpr) == 2 && contains(tmpr{1},'end'))
                        % Subtraction with first term 'end'
                        indr = str2double(tmpr{2});
                        indr = -indr;
                    else
                        error('Group Align: Right Boundary malformed')
                    end
                    inds = [indl indr];
                    % indl and indr might be <=0, then it is counted from end
                    indsnan = isnan(inds);
                    if any(indsnan)
                        error('Group Align input range not recognized: Conversion to numeric value of "%s" failed',...
                            strjoin(margs(indsnan),', '))
                    end
                else
                    error('Group Align: Illegible Range')
                end
                
                C = cell(1,length(l1labels));
                for ls = 1:length(l1labels)
                    strs = strsplit(l1labels{ls},'_');
                    endn = length(strs);
                    inds = (inds <= 0) .* endn + inds;
                    % Convert negative indices, which result from the usage of
                    % 'end' into positive ones
                    curinds = inds(1):inds(2);
                    C{ls} = strjoin(strs(curinds),'_');
                    % Reinjection of all previously removed '_'
                end
                [groupIdentifiers,~,groupalign] = unique(C);
                groupalign = reshape(groupalign,[1 length(groupalign)]);
            else
                [altN,altMatches] = strsplit(groupalign,'alt');
                [groupN,groupMatches] = strsplit(groupalign,'group');
                if ~isempty(altMatches)
                    try 
                        gN = str2double(altN{2});
                        elsperg = length(l1pars)/gN;
                        groupIdentifiers = {groupalign};
                        groupalign = repmat(1:gN,1,elsperg);
                    catch
                        error('Number after "alt" not recognized or not compatible with relto number')
                    end
                elseif ~isempty(groupMatches)
                    try
                        elsperg = str2double(groupN{2});
                        gN = length(l1pars) / elsperg;
                        groupIdentifiers = {groupalign};
                        groupalign = repmat(1:gN,elsperg,1);
                        groupalign = reshape(groupalign,1,length(l1pars));
                    catch
                        error('Number after "group" not recognized or not compatible with relto number')
                    end
                elseif strcmpi(groupalign,'single')
                    groupIdentifiers = {groupalign};
                    groupalign = 1:length(l1pars);
                end
            end
        elseif (checkNumericVector(groupalign) && ...
                length(groupalign) == length(l1pars))
            groupIdentifiers = {'User-Supplied Arrangement'};
        else
            error('Group Align Command not recognized')
        end
        groups = unique(groupalign);
    catch ME
        fprintf('Group Alignment Failed:\n')
        fprintf('%s\n',ME.message)
        return
    end
    
    %% Determining Weights
    
    try
        weights = 1;
        if (isscalar(deform) ...
                && deform ~= 0)
            if isnan(deform)
                metapar = 0.5;
            else
                metapar = deform;
            end
            if (strcmpi(type,'adaptive') ...
                    || strcmpi(type,'AL'))
                
                try
                    arFit(true)
                    weights = 1./abs(ar.p(l1pars)).^metapar;
                catch ME
                    fprintf('OLS Weight assignment not possible, fit failed:\n')
                    fprintf('%s\n',ME.message)
                    return
                end
            elseif strcmpi(type,'adaptive-GW')
                try
                    arFit(true)
                    norms = nan(size(groups));
                    for g = 1:length(groups)
                        ing = l1pars(groupalign == groups(g));
                        norms(g) = norm(ar.p(ing));
                    end
                    weights = 1./norms.^metapar;
                catch ME
                    fprintf('Group-Wise OLS assignment unsuccessful:\n')
                    fprintf('%s\n',ME.message)
                    return
                end
                
            end
        end
    catch ME
        fprintf('Weight assignment not possible:\n')
        fprintf('%s\n',ME.message)
        return
    end
    
    try
        fprintf('Starting Grouped %s-Regularization,\n',type)
        fprintf('Penalizing %i fold-changes,\n',length(l1pars))
        
        tic
        grplasInit(l1pars,means,lbs,ubs,linv,groupalign,weights,...
            p.Results.Threshold,false)
        fprintf('Assembled in %i groups, grouped by\n',ar.grplas.Ngroups)
        fprintf('%s\n',strjoin(groupIdentifiers,', '))
        if strcmpi(lambdascale,'log')
            addstr = '10^';
        else
            addstr = '';
        end
        fprintf('Scanning over %i penalty strengths from %s%.2g to %s%.2g\n',...
            nticks,addstr,lammin,addstr,lammax)
        grplasScan([],[],p.Results.Gradient,[],p.Results.OptimizerSteps)
        grplasUnpen
        grplasSelectOpt([],test)
        if p.Results.doPlot
            grplasPlot
        end
        if p.Results.doTree
            grplasTree
        end
        tmp_scantime = toc;
        
        fprintf('Regularization was successful in %.3f sec.\n',tmp_scantime)
    catch ME
        fprintf('Scan routine for group-wise scan failed:\n')
        fprintf('%s\n',ME.message)
        return
    end
                
end

switch nargout
    case 0
        varargout = {};
    case 1
        varargout = {tmp_scantime};
    case 2
        varargout = {tmp_scantime,find(ar.qFit == 1)};
end


end

