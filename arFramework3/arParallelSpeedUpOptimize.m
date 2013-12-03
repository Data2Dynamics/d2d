function arParallelSpeedUpOptimize(n)

global ar

if(~exist('nruns','var'))
    n = 100;
end
sensis = false;

nReset = ar.config.nParallel;

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
        arChi2(sensis, []);
        arStortTimeConditions(j,0);
    catch err_id
        arStortTimeConditions(j,-1);
        disp(err_id.message);
    end
end
arWaitbar(-1);

labels = cell(1,length(ar.config.threads.ms));
for jt=1:length(ar.config.threads.ms)
    labels{jt} = sprintf('m%i c%i', ar.config.threads.ms(jt)+1, ar.config.threads.cs(jt)+1);
end

arSetParallelThreads(nReset);

%% plot

% sort
% isort = 1:size(ar.config.threads_timings,2);
[~,isort] = sort(median(ar.config.threads_timings), 2, 'descend');
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
title('runtime comparison (sequential)');

ngroups = 2;
S = T(:,1:ngroups);
Sindex = num2cell(1:ngroups);

% recombine
jk = ngroups+1;
while(jk<=size(T,2))
    [~,isort] = sort(median(S));
    S(:,isort(1)) = S(:,isort(1)) + T(:,jk);
    Sindex{isort(1)}(end+1) = jk;
    jk = jk + 1;
end
Sindex{:}

subplot(2,2,2);
boxplot(S, 'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', 'k');
xlabel('computation time (aggregated)');
title('runtime comparison');
xlimtmp = xlim;

X = [];
for j=1:ngroups
    X = [X; ismember(1:size(ar.config.threads_timings,2), Sindex{j})'];
end

%%  optimization

% ga
% p = ones(1,size(T,2));
% Ttmp = ar.config.threads_timings;
% popt = ga(@fkt_ga_opt,length(p),[],[],[],[], ...
%     ones(1,size(Ttmp,2)), ones(1,size(Ttmp,2))*ngroups, [], 1:length(p))
% S2 = fkt_ga_opt2(popt);

nred = 0;

Ttmp = median(ar.config.threads_timings);
btmp = sum(Ttmp) / (ngroups - nred);

A1 = zeros(ngroups,length(Ttmp)*ngroups);
b1 = zeros(1,ngroups);
jk = 1;
for j=1:ngroups
    A1(j,jk:(jk+length(Ttmp)-1)) = Ttmp;
    b1(j) = btmp;
    jk = jk + length(Ttmp);
end

A2 = zeros(ngroups,length(Ttmp)*ngroups);
b2 = zeros(1,ngroups);
for j=1:(length(Ttmp)*ngroups)
    A2(j, mod((1:length(Ttmp)*ngroups)-j,length(Ttmp))==0) = 1;
    b2(j) = 1;
    jk = jk + length(Ttmp);    
end
    
% X = bintprog([],A1,b1);
% X = bintprog([],[],[],A2,b2);
X = bintprog([],A1,b1,A2,b2);

figure(2);
X2 = reshape(X,length(Ttmp),ngroups);
imagesc(X2)

S2 = [];
for j=1:ngroups
    S2(:,j) = sum(ar.config.threads_timings(:,X2(:,j)==1),2);
end

figure(1);
subplot(2,2,4);
boxplot(S2, 'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', 'k');
xlabel('computation time (aggregated)');
title('runtime comparison');
xlim(xlimtmp);

function f = fkt_ga_opt(p)
F = fkt_ga_opt2(p);
f = max(median(F));
% f = median(median(F));

function F = fkt_ga_opt2(p)

global ar

Ttmp = ar.config.threads_timings;
nmax = max(p);

F = zeros(size(Ttmp,1),nmax);
for j = 1:nmax
    F(:,j) = F(:,j) + sum(Ttmp(:,p==j),2);
end


function arStortTimeConditions(j,flag)

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
