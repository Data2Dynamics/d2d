% arIdentifiablityTest(silent, radius, nfit, doFittigFirst)
% 
%   This function implements a fast procedure to check whether structural
%   non-identifiablities are present.
% 
%   This idea is to pull the optimized parameters (penalty for the radius)
%   and check, whether there is any direction where the penalty can be
%   satisfied without deterioate mismatch to the data (e.g. likelihood).
% 
%   THE MODEL HAS TO BE FITTED FIRST
%   Local identifiablity is tested around the optimum.
% 
%   silent   Prevent command-line output ? Default: false
% 
%   radius   target radius when pulling
% 
%   nfit    number of initial guesses/fits (like LHS)
%           random initial guess is drawn in the ball wir uniformly
%           distributed radius [0,radius]
% 
%   doFittingFirst  if true, then the model is first fitted without penalty
%           before arIdentifiablityTest is checked. Could be applied, if it
%           is unknown whether the parameters were already fitted.
%           If true, arFit is called then arFitLHS(nfit-1)
% 
% 
%   Example
% 
% 
%   arLoad(2)
%   arIdentifiablityTest  % default call
% 
%   % userdefined:
%   arIdentifiablityTest([],.1,10)  % 10 random initial guess
% 
% ckreutz 2.6.17

function arIdentifiablityTest(silent, radius, nfit, doFittigFirst)
if ~exist('silent','var') || isempty(silent)
    silent = false;
end
if ~exist('radius','var') || isempty(radius)
    radius = 1;
elseif radius<=0
    error('arIdentifiabilityTest is only reasonable for radius>0')
end
if ~exist('nfit','var') || isempty(nfit)
    nfit = 5;
end
if ~exist('doFittigFirst','var') || isempty(doFittigFirst)
    doFittigFirst = 0;
end



thresh = 1e-3; % threshold for deciding
global ar

if ~silent
    disp('Identifiabilty test started ...');
end

% res_fun = 'user_residual_fun_IdentifiabiltyTest';
res_fun = @user_residual_fun_IdentifiabiltyTest2;

if isfield(ar.config,'user_residual_fun')
    tmp = char(ar.config.user_residual_fun);
    if ~isempty(tmp) && strcmp(tmp,res_fun)~=1
        error('arIdentifiablityTest can not handle user-specific residuals. Remove it via "ar.config.user_residual_fun = '''';"')
    end
end


ar.IdentifiabilityTest = struct;
ar.IdentifiabilityTest.n_rescall = 0;
ar.IdentifiabilityTest.pLabel = ar.pLabel;
ar.IdentifiabilityTest.p0 = ar.p + 0.0;
ar.IdentifiabilityTest.radius = radius;
% ar.IdentifiabilityTest.penaltySD = 1/icdf('chi2',0.95,1);
ar.IdentifiabilityTest.penaltySD = ar.IdentifiabilityTest.radius;
ar.IdentifiabilityTest.thresh  = thresh;
ar.IdentifiabilityTest.nBound0 = sum( (ar.p(ar.qFit==1) - ar.lb(ar.qFit==1)) < ar.config.par_close_to_bound*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)) | ...
    (ar.ub(ar.qFit==1) - ar.p(ar.qFit==1)) < ar.config.par_close_to_bound*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)) ) ;

ar.IdentifiabilityTest.message = cell(0);


%% Berechne objective function without pulling away:
% try
if doFittigFirst
    if ~silent
        disp('Standard fitting ...')
    end
    %     ar.config.optim.Display = 'iter';
    arFit(true)
    if nfit>1
        arFitLHS(nfit-1)
    end
    %     arPrint
end
ar.IdentifiabilityTest.residual_fun = res_fun;
ar.config.user_residual_fun = ar.IdentifiabilityTest.residual_fun;

estimatetime = 0;
tic;

%% without penalty
ar.IdentifiabilityTest.state = 'off';
arCalcMerit(true,ar.p(ar.qFit==1)); %% be sure that the residuals are up-to-date, sensi = true ensures that ODE intergation steps are the same as within fitting
ar.IdentifiabilityTest.chi2_woPenalty0 = arGetMerit(true) + 0.0;


%% before Fitting, with penalty
ar.IdentifiabilityTest.state = 'on';
arCalcMerit(true,ar.p(ar.qFit==1)); %% be sure that the residuals are up-to-date, sensi = true ensures that ODE intergation steps are the same as within fitting
ar.IdentifiabilityTest.chi2_wPenalty0 = arGetMerit(true) + 0.0;

ar.IdentifiabilityTest.message{end+1} = sprintf('\nIdentifiability test was performed with radius = %d and penalty-SD = %d. \n\n',ar.IdentifiabilityTest.radius,ar.IdentifiabilityTest.penaltySD);
ar.IdentifiabilityTest.ps_start = randomParsInBall(nfit,ar.IdentifiabilityTest.radius);
ar.IdentifiabilityTest.ps_start(1,:) = ar.p; % the first one is always ar.p

arFits(ar.IdentifiabilityTest.ps_start);
ar.IdentifiabilityTest.chi2s = ar.chi2s;
ar.IdentifiabilityTest.ps = ar.ps;

ar.IdentifiabilityTest.nBound = sum( (ar.p(ar.qFit==1) - ar.lb(ar.qFit==1)) < ar.config.par_close_to_bound*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)) | ...
    (ar.ub(ar.qFit==1) - ar.p(ar.qFit==1)) < ar.config.par_close_to_bound*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)) ) ;

ar.IdentifiabilityTest.euclDist_relTo_deltaSD = sqrt(sum((ar.p-ar.IdentifiabilityTest.p0).^2))./ar.IdentifiabilityTest.penaltySD;
ar.IdentifiabilityTest.euclDist = sqrt(sum((ar.p-ar.IdentifiabilityTest.p0).^2));
ar.IdentifiabilityTest.p = ar.p + 0.0;
ar.IdentifiabilityTest.dp = ar.IdentifiabilityTest.p - ar.IdentifiabilityTest.p0;
[~,ar.IdentifiabilityTest.dp_order] = sort(-abs(ar.IdentifiabilityTest.dp));


%% after Fitting, with penalty
arCalcMerit(true,ar.p(ar.qFit==1)); %% be sure that the residuals are up-to-date, sensi = true ensures that ODE intergation steps are the same as within fitting
ar.IdentifiabilityTest.chi2_wPenalty = arGetMerit(true) + 0.0;

% ar.IdentifiabilityTest.chi2_ppl = ar.IdentifiabilityTest.chi2_wPenalty - ar.IdentifiabilityTest.penalty;

ar.IdentifiabilityTest.state = 'off';
%% after fitting, without penalty
arCalcMerit(true,ar.p(ar.qFit==1)); %% be sure that the residuals are up-to-date, sensi = true ensures that ODE intergation steps are the same as within fitting
ar.IdentifiabilityTest.chi2_woPenalty = arGetMerit(true) + 0.0;

ar.config.user_residual_fun = '';  % residual_fun für IdentifiablityTest wieder ausschalten
ar.p = ar.IdentifiabilityTest.p0;

% %% calculate a score
% if ar.IdentifiabilityTest.dchi2_woPenalty>=-1e-2
%     ar.IdentifiabilityTest.score = icdf('chi2',ar.IdentifiabilityTest.dchi2_woPenalty,1); % 1 means identifiable, 0 means structurally non-id.
% elseif ar.IdentifiabilityTest.dchi2_woPenalty>-1e-3
%     ar.IdentifiabilityTest.score = 1;
% else
%     ar.IdentifiabilityTest.dchi2_woPenalty
%     warning('-2* Log-likelihood improves in the case of penalized optimization. Possible reasons: Either optimization does not work or the model was not fitted before.')
%     ar.IdentifiabilityTest.score = NaN;
% end




%% Output
ar.IdentifiabilityTest.penalty = ar.IdentifiabilityTest.chi2_wPenalty -ar.IdentifiabilityTest.chi2_woPenalty;
ar.IdentifiabilityTest.penalty0 = ar.IdentifiabilityTest.chi2_wPenalty0 -ar.IdentifiabilityTest.chi2_woPenalty0;
ar.IdentifiabilityTest.dchi2_woPenalty = ar.IdentifiabilityTest.chi2_woPenalty - ar.IdentifiabilityTest.chi2_woPenalty0;
ar.IdentifiabilityTest.dchi2_total = ar.IdentifiabilityTest.chi2_wPenalty - ar.IdentifiabilityTest.chi2_woPenalty0;

ar.IdentifiabilityTest = sortFields(ar.IdentifiabilityTest);

if ar.IdentifiabilityTest.euclDist < 0.01*ar.IdentifiabilityTest.radius
    ar.IdentifiabilityTest.message{end+1} = sprintf('Optimal parameters change less than 1%s of radius. Is optimization working?\n\n','%');
end

if ar.IdentifiabilityTest.nBound0 < ar.IdentifiabilityTest.nBound
    warnmessage = sprintf(' Warning: Penalization force addional parameters to bounds. Decreasing the radius is suggested in this case.\n\n');
    warning(warnmessage)
    ar.IdentifiabilityTest.message{end+1} = warnmessage;
end

chi2sort = sort(ar.IdentifiabilityTest.chi2s);
chi2sort2 = sort(ar.IdentifiabilityTest.chi2s(2:end));
if range(chi2sort)<ar.IdentifiabilityTest.thresh
    ar.IdentifiabilityTest.message{end+1} = sprintf('All %i optimization runs are in the chi2-range %g.\n\n',length(ar.IdentifiabilityTest.chi2s),range(chi2sort));
elseif length(chi2sort2)>2 && range(chi2sort2)<ar.IdentifiabilityTest.thresh
    ar.IdentifiabilityTest.message{end+1} = sprintf('All %i optimization runs with random intial guesses are in the chi2-range %g.\n\n',length(ar.IdentifiabilityTest.chi2s)-1,range(chi2sort2));
elseif sum(chi2sort <(min(chi2sort) + ar.IdentifiabilityTest.thresh)) < length(chi2sort)/2
    anzsame = sum(chi2sort <(min(chi2sort) + ar.IdentifiabilityTest.thresh));
    fracsame = anzsame/length(chi2sort);
    ar.IdentifiabilityTest.message{end+1} = sprintf('Only %i optimization runs (%.2f%s) are in the chi2-range %g. Increasing the number of fits should be considerd. \n\n',anzsame,fracsame*100,'%',ar.IdentifiabilityTest.thresh);
else
    anzsame = sum(chi2sort <(min(chi2sort) + ar.IdentifiabilityTest.thresh));
    fracsame = anzsame/length(chi2sort);
    ar.IdentifiabilityTest.message{end+1} = sprintf('%i optimization runs (%.2f%s) are in the chi2-range %g. \n\n',anzsame,fracsame*100,'%',ar.IdentifiabilityTest.thresh);
end

estimatetime = estimatetime + toc;
ar.IdentifiabilityTest.timing = estimatetime;

ar.IdentifiabilityTest.message{end+1} = sprintf('Calculations took %.2f seconds.\n',ar.IdentifiabilityTest.timing);
if isfield(ar,'ple') && isfield(ar.ple,'timing')
    ar.IdentifiabilityTest.message{end+1} = sprintf('[Compared to %.2f seconds required for calcuating the likelihood profiles.]\n\n',sum(ar.ple.timing));
else
    ar.IdentifiabilityTest.message{end} = sprintf('%s\n',ar.IdentifiabilityTest.message{end});
end

ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (increase of merit by penalty, before fitting)\n',ar.IdentifiabilityTest.chi2_wPenalty0 - ar.IdentifiabilityTest.chi2_woPenalty0);
ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (decrease of merit by fitting)\n',ar.IdentifiabilityTest.chi2_wPenalty0-ar.IdentifiabilityTest.chi2_wPenalty);
ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (movement of parameters by penalized fitting)\n',ar.IdentifiabilityTest.euclDist);
% ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (penalty after fitting)\n',ar.IdentifiabilityTest.penalty);
% ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (increase of data-chi2 by penalty)\n',ar.IdentifiabilityTest.dchi2_woPenalty);
ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (total increase of merit by penalty) PRIMARY CRITERION\n',ar.IdentifiabilityTest.dchi2_total);
% ar.IdentifiabilityTest.message{end+1} = sprintf('%5.4f (movement of parameters rel. to penaltySD)\n',ar.IdentifiabilityTest.euclDist_relTo_deltaSD);


if ar.IdentifiabilityTest.dchi2_total < ar.IdentifiabilityTest.thresh;
    ar.IdentifiabilityTest.message{end+1} = sprintf('\nModel is structurally non-identifiable.\n\n');
    anzNI = sum(abs(ar.IdentifiabilityTest.dp)>0.1*ar.IdentifiabilityTest.radius);
    ar.IdentifiabilityTest.p_nonIdentifiable = ar.pLabel(ar.IdentifiabilityTest.dp_order(1:anzNI));
else
    ar.IdentifiabilityTest.p_nonIdentifiable = cell(0);
    ar.IdentifiabilityTest.message{end+1} = sprintf('\nModel is identifiable.\n\n');
end

if ~silent
    disp(' ')
    for m=1:length(ar.IdentifiabilityTest.message)
        fprintf('%s',ar.IdentifiabilityTest.message{m});
    end
end


function s2 = sortFields(s)
fn = sort(fieldnames(s));
s2 = struct;
for i=1:length(fn)
    s2.(fn{i}) = s.(fn{i});
end




function p0 = randomParsInBall(n,R)
global ar
nqfit = sum(ar.qFit==1);
bolfit = ar.qFit==1;
indfit = find(bolfit);

richtung = rand(n,nqfit);
for i=1:size(richtung,1)
    richtung(i,:) = richtung(i,:) ./ norm(richtung(i,:),2);
end
radius = rand(n,1)*R;

p0 = ones(n,1)*ar.p;
p0(:,bolfit) = p0(:,bolfit) + richtung.*(radius*ones(1,nqfit));

for i=1:size(p0,1)
    bol = p0(i,bolfit)<ar.lb(bolfit);
    p0(i,indfit(bol)) = ar.lb(indfit(bol));
    
    bol = p0(i,bolfit)>ar.ub(bolfit);
    p0(i,indfit(bol)) = ar.ub(indfit(bol));
end


function [res_user, sres_user, res_type] = user_residual_fun_IdentifiabiltyTest2
%   This function penalize the if || ar.p - ar.IdentifiabilityTest.p0 ||_2 
%   unequal to ar.IdentifiabilityTest.radius
%
%   The quadratic difference between || || and radius is penalized, i.e.
%   residuals || || - radius are used by lsqnonlin.

global ar

if strcmp(ar.IdentifiabilityTest.state,'on')
    % fprintf('user_residual_fun_IdentifiabiltyTest2 ...\n');
    res_user = (sqrt(sum((ar.p(ar.qFit==1)-ar.IdentifiabilityTest.p0(ar.qFit==1)).^2)) - ar.IdentifiabilityTest.radius)./ar.IdentifiabilityTest.penaltySD;    
    res_type = 1; % like a data residual
    
    sres_user(1,1:length(ar.p)) = 0;
    if ar.IdentifiabilityTest.radius>0
        zaehler = (2 * ar.p(ar.qFit==1) - 2*ar.IdentifiabilityTest.p0(ar.qFit==1))./ar.IdentifiabilityTest.penaltySD;
        nenner = 2* sqrt(  sum(  (ar.p(ar.qFit==1) - ar.IdentifiabilityTest.p0(ar.qFit==1)).^2)  );
        sres_user(1,ar.qFit==1) = zaehler ./ nenner;
        sres_user(isnan(sres_user)) = 0;
    end
    
    ar.IdentifiabilityTest.n_rescall = ar.IdentifiabilityTest.n_rescall+1;
else % Die Anzahl Datenpkt soll in beiden Fällen gleich sein.
    res_user = 0;
    sres_user = zeros(1,length(ar.p));
    res_type = 1;
end

ar.IdentifiabilityTest.res = res_user;
ar.IdentifiabilityTest.sres = sres_user;
