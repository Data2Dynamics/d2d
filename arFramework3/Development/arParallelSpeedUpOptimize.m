% 
function arParallelSpeedUpOptimize(n)

global ar

if(~exist('nruns','var'))
    n = 100;
end
sensis = false;

if(~isfield(ar.config, 'threads_timings'))
    nReset = ar.config.nThreads;
    
    arSetParallelThreads(1);
    
    rng('shuffle');
    ps = ones(n,1) * ar.p;
    psrand = lhsdesign(n,sum(ar.qFit==1));
    
    psrand = psrand .* (ones(n,1)*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)));
    psrand = psrand + (ones(n,1)*ar.lb(ar.qFit==1));
    
    ps(:,ar.qFit==1) = psrand;
    
    ar.config.threads_timings = nan(n,length(ar.config.threads));
    arWaitbar(0);
    for j=1:n
        arWaitbar(j,n, 'optimizing parallel calculations');
        ar.p = ps(j,:);
        try
            arCalcMerit(sensis, []);
            arSaveTimeConditions(j,0);
        catch err_id
            arSaveTimeConditions(j,-1);
            disp(err_id.message);
        end
    end
    arWaitbar(-1);
    
    labels = cell(1,length(ar.config.threads.ms));
    for jt=1:length(ar.config.threads.ms)
        labels{jt} = sprintf('m%i c%i', ar.config.threads.ms(jt)+1, ar.config.threads.cs(jt)+1);
    end
    
    arSetParallelThreads(nReset);
end

%% plot

nthreads = ar.config.nThreads;

% sort
isort = 1:size(ar.config.threads_timings,2);
% [~,isort] = sort(median_nan(ar.config.threads_timings), 2, 'descend');
if(~exist('labels','var'))
    labels = [];
else
    if(~isempty(labels))
        labels = labels(isort);
    end
end
T = ar.config.threads_timings(:,isort);

figure(1); clf;
subplot(2,2,[1 3]);
boxplot(T, 'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', 'k', 'labels', labels);
xlabel('computation time');
title('runtime sequential');

% standard combination

S0 = zeros(size(T,1),nthreads);

jk = 1;
for jt = 1:size(T,2)
    S0(:,jk) = S0(:,jk) + T(:,jt);
    jk = jk + 1;
    if(jk>nthreads)
        jk = 1;
    end
end

subplot(2,2,2);
boxplot(S0, 'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', 'k');
xlabel('computation time (aggregated)');
title('runtime default');
xlimtmp = xlim;

% recombine

S = T(:,1:nthreads);
Sindex = num2cell(1:nthreads);

jk = nthreads+1;
while(jk<=size(T,2))
    [~,isort] = sort(median_nan(S));
    S(:,isort(1)) = S(:,isort(1)) + T(:,jk);
    Sindex{isort(1)}(end+1) = jk;
    jk = jk + 1;
end
Sindex{:}

subplot(2,2,4);
boxplot(S, 'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', 'k');
xlabel('computation time (aggregated)');
title('runtime optimized');
xlim(xlimtmp);


function y = median_nan(x)

y = nan(1,size(x,2));
for j=1:size(x,2)
    tmp = x(:,j);
    tmp = tmp(~isnan(tmp));
    y(j) = median(tmp);
end

function arSaveTimeConditions(j,flag)

global ar

for jt=1:length(ar.config.threads.ms)
    jm = ar.config.threads.ms(jt) + 1;
    jc = ar.config.threads.cs(jt) + 1;
    
    if(flag==0)
        ar.config.threads_timings(j,jt) = ...
            ar.model(jm).condition(jc).stop - ...
            ar.model(jm).condition(jc).start;
    elseif(flag==1)
        ar.config.threads_timings(j,jt) = nan;
    end
end
