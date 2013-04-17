function arOptimizerTest(range)

if(~exist('range','var'))
    range = -1;
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

N = 30;

mus = linspace(0,10^range,N);

chi2s = {};
chi2s_expect = {};
xs = {};
labels = {'gradient','newton','levenberg-marquardt','trust region'};
method = [3 2 1 0];
styles = {'s-r','o-g','*-b','d-c'};
styles_expect = {'s--r','o--g','*--b','d--c'};

arWaitbar(0);
try %#ok<TRYNC>
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
                ar.p(ar.qFit==1) = p(ar.qFit==1) + dptmp;
                arChi2(true);
                
                chi2s_expect{jm}(j) = llh_expect; %#ok<AGROW>
                xs{jm}(j) = norm(dptmp); %#ok<AGROW>
                chi2s{jm}(j) = ar.chi2fit; %#ok<AGROW>
            end
        end
    end
end
arWaitbar(-1);

ar.p = p;
arChi2(true);

%% plot

figure(1);

for jm = 1:length(method)
    subplot(3,4,jm);
    plot(xs{jm}, chi2s{jm}, styles{jm});
    hold on
    plot(xs{jm}, chi2s_expect{jm}, styles_expect{jm});
    hold off
    title(labels{jm});
end

subplot(3,4,5:12);
for jm = 1:length(method)
    plot(xs{jm}, chi2s{jm}, styles{jm});
    hold on
    plot(xs{jm}, chi2s_expect{jm}, styles_expect{jm});
end
hold off



