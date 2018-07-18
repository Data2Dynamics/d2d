% This function uses the Minimum-Merger Idea for assessing ODE tolerances.
% Bad tolerances yields a jittery curve between two parameter points.
% 
% It turned out, however, that this jitter is NOT a reliable indicator for
% assessing integration tolerances.
% 
% Example:
% tols = logspace(-1,-12,12);
% res = cell(size(tols));
% for i=1:length(tols)
%     ar.config.rtol=tols(i);
%     ar.config.atol=tols(i);
%     arAssessOdeTolerances
%     res{i} = ar.opti_tols;
%     title(['atol = rtol = ',num2str(ar.config.atol)])
% end

function arAssessOdeTolerances

global ar

ar.opti_tols = struct;
ar.opti_tols.np = 100;
ar.opti_tols.dp = 0.01;
ar.opti_tols.maxiter = 20;
ar.opti_tols.ntrial = 10;
ar.opti_tols.chi2s = [];
ar.opti_tols.pstarts = [];
ar.opti_tols.mergers = cell(1,ar.opti_tols.ntrial);

for I=1:ar.opti_tols.ntrial
    %% Fit first with user-defined residual, but state=off:
%     Wenn man arFit nur für I==1 macht und somit ar.p immer auf dem gleichen Wert hat, so kommt immer das Gleiche raus
%     if I==1
        MiniumMergerInit(ar.p, ar.p,   'penaltySD', ar.opti_tols.dp/10, 'maxiter',100); % is only required to have the residual
        ar.opti_tols.merger.state = 'off';
        arFit
        ar.opti_tols.pstart = ar.p;
%     else
%         ar.p = ar.opti_tols.pstart;
%     end
    
    %% Find flat region using the identifiability test:
    arIdentifiablityTest(true,ar.opti_tols.dp);
    ar.opti_tols.pend = ar.opti_tols.pstart + ar.IdentifiabilityTest.dp;
    
    %% Fit along the line:
    ar.opti_tols.pstarts(I,:) = ar.p;
    arMinimumMerger(ar.opti_tols.pstart, ar.opti_tols.pend, 'penaltySD', ar.opti_tols.dp/100, 'maxiter',100);
    ar.opti_tols.chi2s(:,I) = ar.opti_tols.merger.chi2s;
    ar.opti_tols.mergers{I} = ar.opti_tols.merger;
end

% %% Plotting
hold off
plot(ar.opti_tols.chi2s,'.')
xlabel('Fit index 1:N_{radius}')
% % q = quantile(ar.opti_tols.merger.chi2s,[0,.2]);
% % set(gca,'YLim',[q(1)-0.02*diff(q),q(2)+100*diff(q)])





% status = arMinimumMerger(pstart,pend,doplot)
%
%   This functioon merges two minima, i.e. a path from pstart is searched
%   to pend via radial penalization. pend is in the center and the radius
%   is iteratively reduced.

function status = arMinimumMerger(pstart,pend,varargin)
global ar

    MiniumMergerInit(pstart,pend,varargin{:})    
    
    
    pIn = ar.p;
    tolx = ar.config.optim.TolX;
    tolf = ar.config.optim.TolFun;
    maxiter = ar.config.optim.MaxIter;
    
    ar.config.optim.MaxIter = ar.opti_tols.maxiter;    
%     ar.config.optim.TolX = 0;  % always make
    ar.config.optim.TolFun = 0;
    
    try
        
    for i=1:ceil(length(ar.opti_tols.merger.radius)/2)  % for small radius, TolX has to be decreased

        fprintf('.')
        ar.opti_tols.merger.iradius = i;
        ar.opti_tols.merger.ps_start(ar.opti_tols.merger.iradius,:) = ar.p;
        arFit(true)
        ar.opti_tols.merger.ps(ar.opti_tols.merger.iradius,:) = ar.p;
        ar.opti_tols.merger.chi2s(ar.opti_tols.merger.iradius) = arGetMerit(true);
        ar.opti_tols.merger.res_all(ar.opti_tols.merger.iradius) = ar.opti_tols.merger.res;
        
        if i==1
            hold off
        elseif i==2
            hold on
        end
        plot(ar.opti_tols.merger.iradius,ar.opti_tols.merger.chi2s(ar.opti_tols.merger.iradius),'b.');
        if rem(i,5)==0
            drawnow
        end
    end
    fprintf('\n')

    catch ERR
        ar.config.optim.TolX   = tolx;
        ar.config.optim.TolFun = tolf;
        ar.config.optim.MaxIter = maxiter;    

        arRemoveCustomResidual( 'MiniumMerger' ); % remove custom residual
        ar.opti_tols.merger.ps(ar.opti_tols.merger.iradius,:) = ar.p;
        ar.opti_tols.merger.chi2s(ar.opti_tols.merger.iradius) = NaN;
        ar.opti_tols.merger.res_all(ar.opti_tols.merger.iradius) = ar.opti_tols.merger.res;

        ar.p = pIn;
        
        rethrow(ERR)
    end
    
    ar.config.optim.TolX   = tolx;
    ar.config.optim.TolFun = tolf;
    ar.config.optim.MaxIter = maxiter;

    ar.p = pIn;
        
arRemoveCustomResidual( 'MiniumMerger' ); % remove custom residual





% This function is used to initialize merging of minima. The field ar.merge
% is created.
function MiniumMergerInit(pstart,pend,varargin)
if rem(nargin,2) ~= 0
    error('arMiniumMergerInit(pstart,pend, varargin): two arguments have to be provided.')
end

global ar

abst = sqrt(sum((pstart(ar.qFit==1)-pend(ar.qFit==1)).^2));

ar.opti_tols.merger.pstart = pstart;
ar.opti_tols.merger.pend = pend;
ar.opti_tols.merger.resfun = @user_residual_fun_MinimumMerger;
ar.opti_tols.merger.radius = linspace(abst,0,100)';  % 0 läuft noch nicht
% ar.opti_tols.merger.radius = linspace(abst,0,101)';
% ar.opti_tols.merger.radius(end) = []; %remove zero
ar.opti_tols.merger.iradius = 1;  % use the first radius first

ar.opti_tols.merger.maxstep = min(diff(ar.opti_tols.merger.radius))/10;
ar.opti_tols.merger.penaltySD = ar.opti_tols.merger.maxstep;

ar.opti_tols.merger.state = 'on';
ar.opti_tols.merger.n_rescall = 0;
ar.opti_tols.merger.ps_start = NaN(length(ar.opti_tols.merger.radius),length(ar.p));
ar.opti_tols.merger.ps = NaN(length(ar.opti_tols.merger.radius),length(ar.p));
ar.opti_tols.merger.chi2s = NaN(length(ar.opti_tols.merger.radius),1);
ar.opti_tols.merger.maxiter = 20;

for i=2:2:length(varargin)
    ar.opti_tols.merger.(varargin{i-1}) = varargin{i};
end

ar.merger = ar.opti_tols.merger;

arAddCustomResidual( 'MiniumMerger', ar.opti_tols.merger.resfun, 1 );

ar.p = ar.opti_tols.merger.pstart;





function [res_user, res_type, sres_user] = user_residual_fun_MinimumMerger
%   This function penalize the if || ar.p - ar.IdentifiabilityTest.p0 ||_2
%   unequal to ar.IdentifiabilityTest.radius
%
%   The quadratic difference between || || and radius is penalized, i.e.
%   residuals || || - radius are used by lsqnonlin.

global ar

p_target = ar.opti_tols.merger.pend;
radius = ar.opti_tols.merger.radius(ar.opti_tols.merger.iradius);

if strcmp(ar.opti_tols.merger.state,'on')
    % fprintf('user_residual_fun_IdentifiabiltyTest2 ...\n');
    res_user = (sqrt(sum((ar.p(ar.qFit==1)-p_target(ar.qFit==1)).^2)) - radius)./ar.opti_tols.merger.penaltySD;
    res_type = 1; % like a data residual
    
    sres_user(1,1:length(ar.p)) = 0;
    if radius>0
        zaehler = (2 * ar.p(ar.qFit==1) - 2*p_target(ar.qFit==1))./ar.opti_tols.merger.penaltySD;
        nenner = 2* sqrt(  sum(  (ar.p(ar.qFit==1) - p_target(ar.qFit==1)).^2)  );
        sres_user(1,ar.qFit==1) = zaehler ./ nenner;
        sres_user(isnan(sres_user)) = 0;
    end
    
    ar.opti_tols.merger.n_rescall = ar.opti_tols.merger.n_rescall+1;
else % Die Anzahl Datenpkt soll in beiden Fï¿½llen gleich sein.
    res_user = 0;
    sres_user = zeros(1,length(ar.p));
    res_type = 1;
end

ar.opti_tols.merger.res = res_user;
ar.opti_tols.merger.sres = sres_user;







