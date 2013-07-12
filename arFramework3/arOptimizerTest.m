function arOptimizerTest(range, method, N)

if(~exist('range','var'))
    range = -1;
end
if(~exist('method','var'))
    method = [0 6];
end
if(~exist('N','var'))
    N = 10;
end

global ar

p = ar.p + 0;
lb = ar.lb + 0;
ub = ar.ub + 0;

ar.p = p;
arChi2(true);

llh = ar.chi2fit;
sres = ar.sres(:,ar.qFit==1); 
g = -2*ar.res*sres;
H = 2*(sres'*sres);

% H = hessian(@my_hess_fkt,p(ar.qFit==1));
% Ht = load('Ht.mat'); 
% H = Ht.Ht;                                         

fprintf('norm(g) = %f\n', norm(g));
fprintf('cond(H) = %f\n', cond(H));

mus = linspace(0,10^range,N);

chi2s = {};
chi2s_expect = {};
xs = {};
ps = {};

% method:   
%  0 = trust region (based on modified trust.m)
%  1 = Levenberg-Marquardt
%  2 = Newton (with maximal step length mu)
%  3 = gradient descent (with steplength mu)
%  4 = gradient descent (to cauchy point with steplength mu)
%  5 = dogleg
%  6 = generalized trust region (based on modified trust.m)
%  7 = MATLABs trdog
%  8 = Newton pcgr (with maximal step length mu)
%  9 = trdog pcgr (with maximal step length mu)
% 10 = dogleg Newton pcgr
% 11 = dogleg trdog pcgr
% 12 = trdog pcgr (no DM)
% 13 = trdog pcgr (no DG)
% 14 = trdog pcgr Levenberg-Marquardt
% 15 = trdog pcgr 2D subspace 
% 16 = trdog pcgr (no DM) 2D subspace 
labels = arNLSstep;
labels = labels(method+1);

minllhs = [];
arWaitbar(0);
% try %#ok<TRYNC>
    for jm = 1:length(method)
        ar.p = p;
        for j = 1:N
            arWaitbar([jm j],[length(method) N]);
            
            [dptmp, solver_calls, qred, ~, grad_dir_frac, llh_expect] = ...
                arNLSstep(llh, g, H, sres, mus(j), ...
                p(ar.qFit==1), lb(ar.qFit==1), ub(ar.qFit==1), 0, [], 0, method(jm));

            fprintf('%s mu=%8.2g norm(dp)=%8.2g solver_calls=%i grad_dir_frac=%4.2f qred=', ...
                labels{jm}, mus(j), norm(dptmp), solver_calls, grad_dir_frac);
            ired = find(qred);
            for jred = ired
                fprintf('%i ', jred);
            end
            fprintf('\n');
            
            try %#ok<TRYNC>
                ps{jm}(j,:) = dptmp; %#ok<AGROW>
                
                ar.p(ar.qFit==1) = p(ar.qFit==1) + dptmp;
                arChi2(true);
                
                chi2s_expect{jm}(j) = llh_expect; %#ok<AGROW>
                xs{jm}(j) = norm(dptmp); %#ok<AGROW>
                chi2s{jm}(j) = ar.chi2fit; %#ok<AGROW>
            end
        end
        minllhs(jm) = min([chi2s{jm} chi2s_expect{jm}]); %#ok<AGROW>
    end
% end
arWaitbar(-1);

ar.p = p;
arChi2(true);

%% plot

figure(1);
clf
dminllh = llh-minllhs;
dminllh(dminllh<1e-3) = 1e-3;
if(length(method)==1)
    plot(xs{1}, chi2s{1}, 'ko-');
    hold on
    plot(xs{1}, chi2s_expect{1}, 'ko--');
    hold off
    title(labels{1});    
    ylim([llh-dminllh*1.1 llh+dminllh*0.1])
else
    for jm = 1:length(method)
        C1 = arLineMarkersAndColors(jm,[],'o','-');
        C2 = arLineMarkersAndColors(jm,[],'o','--');
        subplot(3,length(method),jm);
        plot(xs{jm}, chi2s{jm}, C1{:});
        hold on
        plot(xs{jm}, chi2s_expect{jm}, C2{:});
        hold off
        title(labels{jm});
        ylim([llh-dminllh(jm)*1.1 llh+dminllh(jm)*0.1])
    end
    
    subplot(3,length(method),(length(method)+1):(3*length(method)));
    for jm = 1:length(method)
        C1 = arLineMarkersAndColors(jm,[],'o','-');
        C2 = arLineMarkersAndColors(jm,[],'o','--');
        plot(xs{jm}, chi2s{jm}, C1{:});
        hold on
        plot(xs{jm}, chi2s_expect{jm}, C2{:});
    end
    plot(xlim, [0 0]+chi2s{1}(1), 'k--');
    hold off
    ylim([llh-max(dminllh)*1.1 llh+max(dminllh)*0.1])
end

% figure(2);
% clf
% for jm = 1:length(method)
%     subplot(1,length(method),jm);
%     h = plot(xs{jm}, ps{jm});
%     title(labels{jm});
%     
%     [~, isort] = sort(std(ps{jm}));
%     isort = isort((end-5):end);
%     h = h(isort);
%     hlabels = {};
%     for j=1:length(h)
%         hlabels{j} = sprintf('%i', isort(j)); %#ok<AGROW>
%     end
%     legend(h,hlabels)
% end

% figure(3);
% subplot(3,1,1)
% plot(g, '*-')
% subplot(3,1,2:3)
% Hmax = max(abs(H(:)));
% imagesc(H, [-Hmax Hmax])
% colorbar


function l = my_hess_fkt(p)
global ar
pRes = ar.p;
try
    ar.p(ar.qFit==1) = p;
    arChi2(true);
    l = ar.chi2fit;
    ar.p = pRes;
catch 
    l = nan;
    ar.p = pRes;
end
