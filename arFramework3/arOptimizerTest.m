function arOptimizerTest(range)

global ar

p = ar.p + 0;

ar.p = p;
arChi2(true);

g = -ar.res*ar.sres;
H = ar.sres'*ar.sres;

% gradient direction
sg = g/norm(g);
distp = -[p-ar.ub; -(p - ar.lb)]; % distance to bounds (should be always positive)
exbounds = [p+sg-ar.ub; -(p+sg - ar.lb)] > 0; % dp bejond bound ?
onbound = distp==0;
qred = sum(onbound & exbounds,1)>0;
if(sum(qred)>0)
    gred = g(~qred);
    sg(:) = 0;
    sg(~qred) = gred/norm(gred);
else
    gred = g;
end

fprintf('norm(g) = %g, reduced by %i\n',norm(gred), sum(qred));

% newton step
sn = transpose(pinv(ar.sres)*ar.res');
% sn = transpose(-inv(H)*g');
% distp = -[p-ar.ub; -(p - ar.lb)]; % distance to bounds (should be always positive)
% exbounds = [p+sn-ar.ub; -(p+sn - ar.lb)] > 0; % dp bejond bound ?
% onbound = distp==0;
% qred = sum(onbound & exbounds,1)>0;
% if(sum(qred)>0)
%     gred = g(~qred);
%     Hred = H(~qred,~qred);
%     sn(:) = 0;
%     sn(~qred) = transpose(-inv(Hred)*gred');
% end

N = 10;

chi2g = nan(1,N);
chi2n = nan(1,N);
chi2t = nan(1,N);

xg = nan(1,N);
xn = nan(1,N);
xt = nan(1,N);

% gradient
xg = linspace(0,10^range,N);
arWaitbar(0);
ar.p = p;
for j = 1:N
    arWaitbar(j,3*N);
    ar.p = p + xg(j)*sg;
    ipp = ar.p<ar.lb | ar.p>ar.ub;
    if(sum(ipp)>0)
        break;
    end
    arChi2(true);
    
    chi2g(j) = ar.chi2fit;
end

% newton
ar.p = p;
xnrel = linspace(0,1,N);
for j = 1:N
    arWaitbar(j+N,3*N);
    dptmp = xnrel(j)*sn;
    ar.p = p + dptmp;
    ipp = ar.p<ar.lb | ar.p>ar.ub;
    if(sum(ipp)>0)
        break;
    end
    try
        arChi2(true);
        
        chi2n(j) = ar.chi2fit;
        xn(j) = norm(dptmp);
    end
end
arWaitbar(-1);

% trust region
ar.p = p;
for j = 1:N
    arWaitbar(j+2*N,3*N);
    dptmp = getStepTrust(g, H, xg(j), p, ar.lb, ar.ub);
    ar.p = p + dptmp;
    ipp = ar.p<ar.lb | ar.p>ar.ub;
    if(sum(ipp)>0)
        break;
    end
    arChi2(true);
    
    chi2t(j) = ar.chi2fit;
    xt(j) = norm(dptmp);
end
arWaitbar(-1);

ar.p = p;
arChi2(true);

figure(1);

subplot(3,3,1);
plot(xg,chi2g,'s-r');
title('gradient');

subplot(3,3,2);
plot(xn,chi2n,'o-g');
title('Newton');

subplot(3,3,3);
plot(xt,chi2t,'*-b');
title('trust region');

subplot(3,3,[4:9]);
plot(xg,chi2g,'s-r');
hold on
plot(xn,chi2n,'o-g');
plot(xt,chi2t,'*-b');
hold off



function dp = getStepTrust(g, H, mu, p, lb, ub)
% generate trial steps
%
% g:    gradient
% H:    Hessian matrix
% mu:   trust region size
% p:    current parameters
% lb:   lower bounds
% ub:   upper bounds

dp = trust(-g',H,mu)'; % trust region solution

% PROBLEM: ensuring norm(dp)<=mu in trust function
if(norm(dp)>mu)
    dp = dp/norm(dp)*mu;
end

distp = -[p-ub; -(p - lb)]; % distance to bounds (should be always positive)
exbounds = [p+dp-ub; -(p+dp - lb)] > 0; % dp bejond bound ?

% if p on bounds and dp pointing bejond bound, reduce problem
onbound = distp==0;
qred = sum(onbound & exbounds,1)>0;
if(sum(~qred)==0)
    error('solution outside bounds, no further step possible');
end
if(sum(qred)>0)
    gred = g(~qred);
    Hred = H(~qred,~qred);
    dp_red = trust(-gred',Hred,mu)'; % reduced trust region solution
    
    % PROBLEM: ensuring norm(dp)<=mu in trust function
    if(norm(dp_red)>mu)
        dp_red = dp_red/norm(dp_red)*mu;
    end
    
    dp(:) = 0;
    dp(~qred) = dp_red;
    exbounds = [p+dp-ub; -(p+dp - lb)] > 0; % dp bejond bound ?
end

% if dp too long cut to bounds
dptmp = [dp; dp];
dpredfac = abs(min(distp(exbounds)./dptmp(exbounds)));
if(~isempty(dpredfac))
    if(dpredfac==0)
        error('zero step size');
    end
    dp = dp * dpredfac;
end
