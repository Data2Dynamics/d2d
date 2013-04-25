function arOptimizerTest(range, imethod, N)

if(~exist('range','var'))
    range = -1;
end
if(~exist('imethod','var'))
    imethod = 1:5;
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
g = -2*ar.res*ar.sres(:,ar.qFit==1);
H = 2*ar.sres(:,ar.qFit==1)'*ar.sres(:,ar.qFit==1);

fprintf('norm(g) = %f\n', norm(g));
fprintf('cond(H) = %f\n', cond(H));

mus = linspace(0,10^range,N);

chi2s = {};
chi2s_expect = {};
xs = {};
ps = {};
labels = {'trust region','levenberg-marquardt','newton','gradient','dogleg'};
method = [0 1 2 3 4];
styles = {'d-c','*-b','o-g','s-r','x-m'};
styles_expect = {'d--c','*--b','o--g','s--r','x--m'};

labels = labels(imethod);
method = method(imethod);
styles = styles(imethod);
styles_expect = styles_expect(imethod);

arWaitbar(0);
% try %#ok<TRYNC>
    for jm = 1:length(method)
        ar.p = p;
        for j = 1:N
            arWaitbar([jm j],[length(method) N]);
            [dptmp, solver_calls, qred, ~, grad_dir_frac, llh_expect] = arNLSstep(llh, g, H, mus(j), p(ar.qFit==1), lb(ar.qFit==1), ub(ar.qFit==1), ...
                0, [], 0, method(jm));
            fprintf('%s mu=%4.2f solver_calls=%i grad_dir_frac=%4.2f qred=', ...
                labels{jm}, mus(j), solver_calls, grad_dir_frac);
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
    end
% end
arWaitbar(-1);

ar.p = p;
arChi2(true);

%% plot

figure(1);
clf

if(length(method)==1)
    plot(xs{1}, chi2s{1}, styles{1});
    hold on
    plot(xs{1}, chi2s_expect{1}, styles_expect{1});
    hold off
    title(labels{1});
else
    for jm = 1:length(method)
        subplot(3,length(method),jm);
        plot(xs{jm}, chi2s{jm}, styles{jm});
        hold on
        plot(xs{jm}, chi2s_expect{jm}, styles_expect{jm});
        hold off
        title(labels{jm});
    end
    
    subplot(3,length(method),(length(method)+1):(3*length(method)));
    for jm = 1:length(method)
        plot(xs{jm}, chi2s{jm}, styles{jm});
        hold on
        plot(xs{jm}, chi2s_expect{jm}, styles_expect{jm});
    end
    plot(xlim, [0 0]+chi2s{1}(1), 'k--');
    hold off
end

figure(2);
clf
for jm = 1:length(method)
    subplot(1,length(method),jm);
    h = plot(ps{jm});
    title(labels{jm});
    
    [~, isort] = sort(std(ps{jm}));
    isort = isort((end-5):end);
    h = h(isort);
    hlabels = {};
    for j=1:length(h)
        hlabels{j} = sprintf('%i', isort(j)); %#ok<AGROW>
    end
    legend(h,hlabels)
end

figure(3);
subplot(3,1,1)
plot(g, '*-')
subplot(3,1,2:3)
Hmax = max(abs(H(:)));
imagesc(H, [-Hmax Hmax])
colorbar

