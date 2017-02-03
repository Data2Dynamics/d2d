% perform MOTA on Profile Likelihood Exploit
% to reveal functionally related parameters
%
% pleMOTA

function pleMOTA

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end

% Setup k-Matrix
ps = [];
for j=1:length(ar.ple.ps)
    if(~isempty(ar.ple.ps{j}))
        q = ~isnan(ar.ple.chi2s{j});
        ps = [ps;ar.ple.ps{j}(q,ar.qFit)];
    end
end
    
%% call MOTA
T1 = 0.01; %default = 0.01
T2 = 0.07; %default = 0.07
T3 = 0.08; %default = 0.08

if(size(ps,1) > 300)
    sampleSize = 100;
else
    sampleSize = floor(size(ps,1)/2);
end
    
numOfBootSamp = 30; % default = 35
responseCut = length(ps(1,:)) - 1;

[S, rSquared] = mota(ps, 0, T1, T2, T3, sampleSize, numOfBootSamp, responseCut);
h = motaR(ps, S, rSquared);
ar.ple.mota_results = makeTable(h);

print(h)
% r2 : multiple-r^2 / coefficient of determination
%       ss_res=sum((RESPONSE-sum(PREDICTOR,2)).^2) (sum of squared residuals)
%       ss_tot=sum((RESPONSE-mean(RESPONSE)).^2);
%       r2=1-(ss_res/ss_tot);
%       linear regression: r = correlation
% WIKI: r^2 is a statistic that will give some information about the 
%       goodness of fit of a model.
%       (http://en.wikipedia.org/wiki/Coefficient_of_determination)
%
% cv : std(ps) ./ abs(mean(ps)) %% MEDIAN would be better?
%
% #  : number of same relations found

count = 1;
for j=1:length(ar.ple.p)
    if(ar.qFit(j))
        fprintf('p%i : %s\n', count, ar.ple.p_labels{j});
        count = count + 1;
    end
end

%% Plot

% set relation matrix
ar.ple.IDrelations = nan(length(ar.ple.p));
ar.ple.IDrelations(ar.qFit,ar.qFit) = ...
    ar.ple.mota_results(1:sum(ar.qFit),1:sum(ar.qFit));

% cv influence
ar.ple.IDcv = nan(1,length(ar.ple.p));
ar.ple.IDcv(ar.qFit) = ar.ple.mota_results(:,sum(ar.qFit)+2);
% ar.ple.IDcv(ar.qFit) = abs((quantile(ps,0.75) - quantile(ps,0.15))./quantile(ps,0.5));

% r^2
ar.ple.IDr2 = nan(1,length(ar.ple.p));
ar.ple.IDr2(ar.qFit) = ar.ple.mota_results(:,sum(ar.qFit)+1);

IDweigthing = ar.ple.IDr2; 
% IDweigthing = ar.ple.IDcv;

tmprel = zeros([size(ar.ple.IDrelations) 3]);
for j=1:length(ar.ple.IDstatus)
    if(ar.qFit(j))
        if(ar.ple.IDstatus(j) == 1)
            tmprel(j,:,2) = ar.ple.IDrelations(j,:).*IDweigthing(j);
        elseif(ar.ple.IDstatus(j) == 2)
            tmprel(j,:,2) = ar.ple.IDrelations(j,:).*IDweigthing(j);
            tmprel(j,:,1) = ar.ple.IDrelations(j,:).*IDweigthing(j);
        elseif(ar.ple.IDstatus(j) == 3)
            tmprel(j,:,1) = ar.ple.IDrelations(j,:).*IDweigthing(j);
        elseif(ar.ple.IDstatus(j) == 4)
            tmprel(j,:,3) = ar.ple.IDrelations(j,:).*IDweigthing(j);
        end
        if(~isnan(ar.ple.IDstatus(j)))
            tmprel(j,j,tmprel(j,j,:)~=0) = 1;
        else
            tmprel(j,j,:) = 1;
        end
    end
end
tmprel(tmprel>1) = 1;
tmprel(~ar.qFit, :, :) = 0.5;
tmprel(:, ~ar.qFit, :) = 0.5;

figure(figure(length(ar.ple.p)+1))
image(tmprel);
set(gca, 'YTick', 1:length(ar.ple.p))
set(gca, 'YTickLabel', ar.ple.p_labels)
set(gca, 'XTick', 1:length(ar.ple.p))
set(gca, 'XTickLabel', ar.ple.p_labels)
%rotateticklabel(gca,45)
title('MOTA parameter relations')

% save
saveas(gcf, [ar.ple.savePath '/mota_relations'], 'fig')
saveas(gcf, [ar.ple.savePath '/mota_relations'], 'png')

%% save
pleSave(ar.ple)

