% This function tries optimization step caluclation strategies using
% several submethods of arNLSstep 
% 
function arOptimizerTest(range, submethod, N)

if(~exist('range','var'))
    range = -1;
end
if(~exist('submethod','var'))
    submethod = [0];
end
if(~exist('N','var'))
    N = 10;
end

global ar

p = ar.p + 0;
lb = ar.lb + 0;
ub = ar.ub + 0;

ar.p = p;
arCalcMerit(true,[]);

llh = arGetMerit('chi2fit') + arGetMerit('chi2constr');
res = [ar.res ar.constr];
if(~isempty(ar.sconstr))
    sres = [ar.sres(:,ar.qFit==1); ar.sconstr(:,ar.qFit==1)];
else
    sres = ar.sres(:,ar.qFit==1);
end

g = -2*res*sres;
H = 2*(sres'*sres);

% H = hessian(@my_hess_fkt,p(ar.qFit==1));
% Ht = load('Ht.mat'); 
% H = Ht.Ht;                                         

fprintf('norm(g) = %g\n', norm(g));
fprintf('cond(H) = %g\n', cond(H));

% mus = linspace(0,10^range,N);
mus = [0 logspace(range-5,range,N-1)];

chi2s = {};
chi2s_expect = {};
xs = {};
ps = {};

% submethod:   
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
labels = labels(submethod+1);

minllhs = [];
arWaitbar(0);
% try %#ok<TRYNC>
    for jm = 1:length(submethod)
        ar.p = p;
        for j = 1:N
            arWaitbar([jm j],[length(submethod) N]);
            
            [dptmp, solver_calls, qred, grad_dir_frac, llh_expect] = ...
                arNLSstep(llh, g, H, sres, mus(j), ...
                p(ar.qFit==1), lb(ar.qFit==1), ub(ar.qFit==1), 0, [], submethod(jm));
            
            fprintf('%s mu=%8.2g norm(dp)=%8.2g solver_calls=%i grad_dir_frac=%4.2f cond(H)=%g qred=', ...
                labels{jm}, mus(j), norm(dptmp), solver_calls, grad_dir_frac, cond(H(qred==0,qred==0)));
            ired = find(qred);
            for jred = ired
                fprintf('%i ', jred);
            end
            fprintf('\n');
            
            try 
                ps{jm}(j,:) = dptmp; %#ok<AGROW>
                
                ar.p(ar.qFit==1) = p(ar.qFit==1) + dptmp;
                arCalcMerit(true,[]);
                
                chi2s_expect{jm}(j) = llh_expect; %#ok<AGROW>
                xs{jm}(j) = norm(dptmp); %#ok<AGROW>
                chi2s{jm}(j) = arGetMerit('chi2fit') + arGetMerit('chi2constr'); %#ok<AGROW>
            catch err_id
                disp(err_id.message)
                chi2s_expect{jm}(j) = nan; %#ok<AGROW>
                xs{jm}(j) = nan; %#ok<AGROW>
                chi2s{jm}(j) = nan; %#ok<AGROW>
            end
        end
        minllhs(jm) = min([chi2s{jm} chi2s_expect{jm}]); %#ok<AGROW>
    end
% end
arWaitbar(-1);

ar.p = p;
arCalcMerit(true,[]);

%% plot

figure(1);
clf
dminllh = llh-minllhs;
dminllh(dminllh<1e-3) = 1e-3;
if(length(submethod)==1)
    plot(xs{1}, chi2s{1}, 'ko-');
    hold on
    plot(xs{1}, chi2s_expect{1}, 'ko--');
    hold off
    title(labels{1});    
    ylim([llh-dminllh*1.1 llh+dminllh*0.1])
else
    for jm = 1:length(submethod)
        C1 = arLineMarkersAndColors(jm,length(submethod),[],'o','-');
        C2 = arLineMarkersAndColors(jm,length(submethod),[],'o','--');
        subplot(3,length(submethod),jm);
        plot(xs{jm}, chi2s{jm}, C1{:});
        hold on
        plot(xs{jm}, chi2s_expect{jm}, C2{:});
        hold off
        title(labels{jm});
        ylim([llh-dminllh(jm)*1.1 llh+dminllh(jm)*0.1])
    end
    
    subplot(3,length(submethod),(length(submethod)+1):(3*length(submethod)));
    for jm = 1:length(submethod)
        C1 = arLineMarkersAndColors(jm,length(submethod),[],'o','-');
        C2 = arLineMarkersAndColors(jm,length(submethod),[],'o','--');
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
% for jm = 1:length(submethod)
%     subplot(1,length(submethod),jm);
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
    arCalcMerit(true,[]);
    l = arGetMerit('chi2fit');
    ar.p = pRes;
catch 
    l = nan;
    ar.p = pRes;
end
