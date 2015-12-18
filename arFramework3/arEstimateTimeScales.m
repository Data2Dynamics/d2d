% TS = arEstimateTimeScales(opt,doplot)
% 
%   Estimates the (maximal, i.e. slowest) time scale of the dynamics using
%   the transient function.
%   The time scale is the decrease to 1/e.
% 
% upper bound for TS ist set to ar.model(m).tLim(2)*10
% 
% Examples:
% 
% % 1) Time scale of the data:
% TS = arEstimateTimeScales
% 
% % 2) Time scale of the model simulation (be sure that the simu was updated
% % for the current parameters!)
% TS = arEstimateTimeScales('x');
% 
% 3) incl. plotting:
% TS = arEstimateTimeScales([],true);

function TS = arEstimateTimeScales(opt,doplot)
global ar
if(~exist('opt','var') || isempty(opt))
    opt = '';
end
if(~exist('doplot','var'))
    doplot = false;
end

switch opt
    case 'x'
        fn = 'yExpSimu';
        fn2 = 'ystdExpSimu';
    otherwise  % use data
        fn = 'yExp';
        if ar.config.fiterrors ==1
            fn2 = 'ystdExpSimu';
        else
            fn2 = 'yExpStd';
        end
end

TS = cell(size(ar.model));
for m=1:length(ar.model)
    if(isfield(ar.model(m),'yNames'));
        ynames = ar.model(m).yNames;
    else
        ynames = unique([ar.model(m).data.yNames]);
    end
    TS{m} = NaN(length(ar.model(m).data),length(ynames));
    
    for d=1:length(ar.model(m).data)
        y = ar.model(m).data(d).(fn);
        t = ar.model(m).data(d).tExp;
        sig = ar.model(m).data(d).(fn2);
    
        doy = find(sum(~isnan(y) & abs(y)>10*eps,1)>4);        
        for i=1:length(doy)  
            try
                erg = fitting(t,y(:,doy(i)),sig(:,doy(i)));
                iy = strmatch(ar.model(m).data(d).yNames{doy(i)},ynames,'exact');
                % Model: f(t) = p4*(1-exp(-p1*t)).*exp(-p2*t) + p3*(1-exp(-p1*t))
                % pfit contains est. params at unlog-scale.
                    
                if(erg.pfit(4)>erg.pfit(3)) % there is a second, transient time scale
%                   % the following was tested for a=0.1, a=0.2 but a=1 was best.
                    % if(erg.pfit(4)>a*erg.pfit(3))
                    
%                     TS{m}(d,iy)= min(1/erg.pfit(2),ar.model(m).tLim(2)*10);
                    % some examples suggest, that the sum of both times are
                    % close to what we want:
                    TS{m}(d,iy)= min(sum(1./erg.pfit(1:2)),ar.model(m).tLim(2)*10); % only a sustained signal
                else
                    TS{m}(d,iy)= min(1/erg.pfit(1),ar.model(m).tLim(2)*10); % only a sustained signal
                end
                
                if doplot
                    plotting(erg)
                    hold on
                    plot(TS{m}(d,iy)*ones(1,2),ylim,'k--','LineWidth',2)
                    xl = xlim;
                    xlim([0,xl(2)])
                    input('press a button')
                end
            catch
                disp(lasterr)
            end
        end
    end
end


end




function erg = fitting(t,y,sig)
doneg = true; %logical specifying whether transient deactivation should also be fitted

fun = @model_SustainedVsTransient_WithOffset;
tlev = unique(t);
p0 = [0.5/(tlev(2)-tlev(1)),.5/range(t)/2,.7,0.7,0];
dolog = 1:4;
dp = 1e3;

if(~exist('sig','var') | isempty(sig))
    sig = ones(size(t))*range(y(:))*.1;
end

if(size(t,2)==1)
    t = t'; % make row
end
if(size(y,2)==1)
    y = y'; % make row
end
if(size(sig,2)==1)
    sig = sig'; % make row
end
if(size(p0,2)==1)
    p0 = p0'; % make row
end

if(size(y,1)>1 && size(t,1)==1)
    t = ones(size(y,1),1)*t;  % duplicate rows
end
if(size(y,1)>1 && size(sig,1)==1)
    sig = ones(size(y,1),1)*sig;  % duplicate rows
end
if(size(y,1)>1 && size(p0,1)==1)
    p0 = ones(size(y,1),1)*p0;  % duplicate rows
end

normfak = range(quantile(y,[.2,.8]));
ynorm = y./normfak;


erg.y = y;
erg.t = t;
erg.fun = fun;
erg.sig = sig;
erg.dolog = dolog;
erg.p0 = p0;

erg.pest = NaN(size(erg.y,1),size(p0,2));
erg.pfit = NaN(size(erg.y,1),size(p0,2));
erg.lb = NaN(size(erg.y,1),size(p0,2));
erg.ub = NaN(size(erg.y,1),size(p0,2));
erg.fest = NaN(size(erg.y,1),1);
erg.signum = NaN(size(erg.y,1),1);

erg.tplot = NaN(size(erg.y,1),201);
erg.yplot = NaN(size(erg.y,1),201);
leer = [];
for i=1:size(y,1)    
    % positive direction
    bound = feval(fun,t(i,:),leer,ynorm(i,:)); % boundaries of reasonable parameters
    
    lb = p0(i,:)-dp;
    lb(dolog) = p0(i,dolog)/dp;
    lb = max([lb;bound.lb]);
    ub = p0(i,:)+dp;
    ub(dolog) = p0(i,dolog)*dp;
    ub = min([ub;bound.ub]);

    p0_log = p0(i,:);
    p0_log(dolog) = log(p0(i,dolog));

    lb_log = lb;
    lb_log(dolog) = log(lb(dolog));
    ub_log = ub;
    ub_log(dolog) = log(ub(dolog));

    out = find(p0_log > ub_log | p0_log < lb_log);
    p0_log(out) = mean([lb_log(out);ub_log(out)],1);

    indnotnan = find(~isnan(y(i,:)));
   [ppos,fpos,res,flag,out] = lsqnonlin('assess_lsqnonlin_ln_dolog',p0_log,lb_log,ub_log,...
        optimset('Display','off','Jacobian','on','DerivativeCheck','off'),fun,t(i,indnotnan),ynorm(i,indnotnan),sig(i,indnotnan),dolog);

    if(doneg)
        % negative direction
        boundNeg = feval(fun,t(i,:),leer,-ynorm(i,:)); % boundaries of reasonable parameters

        lbNeg = p0(i,:)-dp;
        lbNeg(dolog) = p0(i,dolog)/dp;
        lbNeg = max([lbNeg;boundNeg.lb]);
        ubNeg = p0(i,:)+dp;
        ubNeg(dolog) = p0(i,dolog)*dp;
        ubNeg = min([ubNeg;boundNeg.ub]);

        p0_logNeg = p0(i,:);
        p0_logNeg(dolog) = log(p0(i,dolog));

        lb_logNeg = lbNeg;
        lb_logNeg(dolog) = log(lbNeg(dolog));
        ub_logNeg = ubNeg;
        ub_logNeg(dolog) = log(ubNeg(dolog));

        out = find(p0_logNeg > ub_logNeg | p0_logNeg < lb_logNeg);
        p0_logNeg(out) = mean([lb_logNeg(out);ub_logNeg(out)],1);

        [pneg,fneg,res,flag,out] = lsqnonlin('assess_lsqnonlin_ln_dolog',p0_logNeg,lb_logNeg,ub_logNeg,...
            optimset('Display','off','Jacobian','on','DerivativeCheck','off'),fun,t(i,indnotnan),-ynorm(i,indnotnan),sig(i,indnotnan),dolog);
    else
        fneg = Inf;
    end
    
    if(fpos < fneg) % choose the better fit
        erg.lb(i,:) = lb;
        erg.ub(i,:) = ub;
        erg.fest(i) = fpos;

        erg.pest(i,:)= ppos;
        erg.pfit(i,:) = erg.pest(i,:);
        erg.pfit(i,dolog) = exp(erg.pest(i,dolog));

        erg.signum(i) = 1;
        erg.yfit(i,:) = feval(fun,t(i,:),erg.pfit(i,:))';
    else
        erg.lb(i,:) = lbNeg;
        erg.ub(i,:) = ubNeg;
        erg.fest(i) = fneg;

        erg.pest(i,:)= pneg;
        erg.pfit(i,:) = erg.pest(i,:);
        erg.pfit(i,dolog) = exp(erg.pest(i,dolog));
        
        erg.signum(i) = -1;
        erg.yfit(i,:) = -feval(fun,t(i,:),erg.pfit(i,:))';
    end

    erg.pfit(i,3:5) = erg.pfit(i,3:5)*normfak;
    
    erg.tplot(i,:) = linspace(nanmin(erg.t(i,:)),nanmax(erg.t(i,:)),201);
    erg.yplot(i,:) = erg.signum(i)* feval(erg.fun,erg.tplot(i,:),erg.pfit(i,:))';
end

end



% Sustainter Anteil hat die schnelle Zeitskala
% 
% f(t) = p4*(1-exp(-p1*t)).*exp(-p2*t) + p3*(1-exp(-p1*t));
% 
% p(1) = fast time-scale, occurs in sustained and transient part
% p(2) = slow transient time-scale
% p(3) = magnitude sustained
% p(4) = magnitude transient
% p(5) = offset

function [y,dydp] = model_SustainedVsTransient_WithOffset(t,p,data)
if(isempty(p))  % return constraints for the parameters
    t = unique(t);
    D = range(data);
    y.lb = [.1/(t(end)-t(1)),.01/(t(end)-t(1)),1e-10,1e-10,  nanmin(data)-0.5*D];
    y.ub = [2/(t(2)-t(1)),1/(t(2)-t(1)),D*2,D*2,  nanmax(data)];
else
    if(size(t,1)==1)
        t = t';
    end

    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    p4 = p(4);

    y1 = p4*(1-exp(-p1*t)).*exp(-p2*t);
    y2 = p3*(1-exp(-p1*t));
    y = y1+y2+p(5);

    if(nargout>1)
        dydp1 = p4*t.*exp(-p1*t).*exp(-p2*t)+p3*t.*exp(-p1*t);
        dydp2 = -p4*(1-exp(-p1*t)).*t.*exp(-p2*t);
        dydp3 = 1-exp(-p1*t);
        dydp4 = (1-exp(-p1*t)).*exp(-p2*t);
        dydp5 = ones(size(t));

        dydp = [dydp1,dydp2,dydp3,dydp4,dydp5];
    end
end

end


% Calculation of the maxumim of the transient part:
% F = 'p4*(1-exp(-p1*t))*exp(-p2*t)';
% dfdt = diff(F,'t')
% tmax = char(solve([char(dfdt),'=0'],'t'))
% Fmax = char(subs(F,'t',tmax))
%   p4*(1-exp(-p1*(-log(p2/(p1+p2))/p1)))*exp(-p2*(-log(p2/(p1+p2))/p1))
% 
% Integrals:
% int('p4*(1-exp(-p1*t))*exp(-p2*t)','t')
% int('p3*(1-exp(-p1*t))','t')

function ass = assess_transient(p,t)
if(sum(p(1:4)<0)>0)
    disp('assess_transient_f1.m, p negative, log considered?')
end
p1 = p(1);
p2 = p(2);
p3 = p(3);
p4 = p(4);

trans = p4*(1-exp(-p1*(-log(p2/(p1+p2))/p1)))*exp(-p2*(-log(p2/(p1+p2))/p1));
sus = p3;

ass.sus = sus/(trans+sus);
ass.trans = trans/(trans+sus);

tend = nanmax(t);
t0 = nanmin(t);

% susInt = p3*(tend+1/p2*exp(-p2*tend)) ...
%     - p3*(t0+1/p2*exp(-p2*t0));

susInt = p3*(1/2*tend^2-1/p1^2*exp(-p1*tend)) - p3*(1/2*t0^2-1/p1^2*exp(-p1*t0));

transInt = p4*(-1/p2/exp(p2*tend)-1/(-p2-p1)*exp(-p2*tend-p1*tend))...
    -p4*(-1/p2/exp(p2*t0)-1/(-p2-p1)*exp(-p2*t0-p1*t0));

ass.susInt = susInt/(susInt+2*transInt);
ass.transInt = 2*transInt/(susInt+2*transInt);

end

function plotting(erg,iplot)
if(~exist('iplot','var') | isempty(iplot))
    iplot = 1:min(10,size(erg.pfit,1)); % the first 10 plots
end
close all

for i=iplot %:min(maxplots,size(erg.pfit,1))
    figure
    plot(erg.tplot(i,:),erg.yplot(i,:),'-');
    hold on    
    errorbar(erg.t(i,:),erg.y(i,:),erg.sig(i,:).*ones(size(erg.y(i,:))),'o')
    
    if(isfield(erg,'ass'))
        title(sprintf('Magnitude: %.2f%s sustained, integral: %.2f%s sustained',100*erg.ass(i).sus,'%',100*erg.ass(i).susInt,'%'))
    end
end

end


% [res,J] = assess_leastsquares_ln(p,model_f,t,y,sigma,varargin)
% 
%  Least-Squares
% 
%   p           parameter array (which is optimized)
%   model_f     the model leading to the model predictions by applying
%                      feval(model_f,t,p,varargin{:})
%   t           times, predictor variables passed to the model
%   y           data, used (only) in this function
%   sigma       amount of noise (SD), used (only) in this function
%   varargin    further arguments passed to the model


function [res,J,H] = assess_lsqnonlin_ln_dolog(p,model_f,t,y,sigma,dolog,varargin)
if(size(t,1)==1)
    t = t';
end
if(size(y,1)==1)
    y = y';
end
if(~exist('sigma','var') | isempty(sigma))
    sigma = ones(size(y));
end
if(~exist('dolog','var') | isempty(dolog))
    dolog = 1:length(p);
end
if(size(sigma,1)==1)
    sigma = sigma';
end

if(size(p,2)==1)
    p = p';
end
punlog = p;
punlog(dolog) = exp(p(dolog));


if(nargout>1)
    dpun_dp = ones(size(p));
    dpun_dp(dolog) = exp(p(dolog));
    
    [yobs,dyobs_dp] = feval(model_f,t,punlog,varargin{:});
    dyobs_dp = dyobs_dp.*(ones(size(dyobs_dp,1),1)*dpun_dp);

    if(max(size(sigma))==1)
        sigma = sigma*ones(size(dyobs_dp,1),1);
    end

    Sens = dyobs_dp./(sigma*ones(1,size(dyobs_dp,2)));
    if(nargout>2)
        H = 2*Sens' * Sens;  % ndop x ndop
    end

    J = 2*dyobs_dp./(sigma*ones(1,length(p)));
else
    yobs = feval(model_f,t,punlog,varargin{:});
end
res = 2*(yobs-y)./(sigma);  % lsqnonlin optimiert 1/2* (res).^2

end