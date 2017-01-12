function arOptimizerTestApprox(range, N, ips, H)

if(~exist('range','var'))
    range = -1;
end
if(~exist('N','var'))
    N = 10;
end
if(~exist('ips','var'))
    ips = 1;
end

global ar

p = ar.p + 0;

arCalcMerit(true,[]);

llh = ar.chi2fit + ar.chi2constr;
res = [ar.res ar.constr];
sres = [ar.sres; ar.sconstr]; 

g = -2*res*sres;
if(~exist('H','var'))
    H = 2*(sres'*sres);
end

fprintf('norm(g) = %g\n', norm(g(ar.qFit==1)));
fprintf('cond(H) = %g\n', cond(H(ar.qFit==1,ar.qFit==1)));

abs(g)

mus = [-fliplr(linspace(0,10^range,N-1)) 0 linspace(0,10^range,N-1)];
% mus = [-fliplr(logspace(range-5,range,N-1)) 0 logspace(range-5,range,N-1)];

chi2s = {};
chi2s_expect = {};

arWaitbar(0);
for jm = 1:length(ips)
    for j = 1:length(mus)
        arWaitbar([jm j],[length(ips) length(mus)]);
        
        try
            dp = zeros(size(p));
            dp(ips(jm)) = mus(j);
            
            llh_expect = llh - g*dp' + 0.5*dp*H*dp';
            chi2s_expect{ips(jm)}(j) = llh_expect; %#ok<AGROW>
            
            ar.p = p + dp;
            arCalcMerit(true,[]);
            chi2s{ips(jm)}(j) = ar.chi2fit + ar.chi2constr; %#ok<AGROW>
        catch err_id
            disp(err_id.message)
            chi2s{ips(jm)}(j) = nan; %#ok<AGROW>
        end
    end
end
arWaitbar(-1);

ar.p = p;
arCalcMerit(true,[]);

%% plot

figure(1)

h = [];
lab = {};
for jm = 1:length(ips)
    C1 = arLineMarkersAndColors(jm,length(ips),[],'o','-');
    C2 = arLineMarkersAndColors(jm,length(ips),[],'o','--');
        
    subplot(1,2,1)
    plot(mus, chi2s{ips(jm)}, C1{:});
    hold on
    plot(mus, chi2s_expect{ips(jm)}, C2{:});
    
    subplot(1,2,2)
    h(jm) = plot(mus, chi2s{ips(jm)}./chi2s_expect{ips(jm)}, C1{:});
    lab{jm} = strrep(ar.pLabel{ips(jm)}, '_' ,'\_');
    hold on
end
subplot(1,2,1)
plot([0 0], ylim, 'k--');
hold off
subplot(1,2,2)
plot(xlim, [1 1], 'k--');
plot([0 0], ylim, 'k--');
hold off
legend(j, lab);

