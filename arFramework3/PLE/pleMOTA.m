% perform MOTA on Profile Likelihood Exploit
% to reveal functionally related parameters
%
% pleMOTA

function pleMOTA

global pleGlobals;

if(isempty(pleGlobals))
    error('perform ple before usage');
end
if(isempty(pleGlobals.ps))
    return
end

% Setup k-Matrix
ps = [];
for j=1:length(pleGlobals.ps)
    if(~isempty(pleGlobals.ps{j}))
        q = ~isnan(pleGlobals.chi2s{j});
        ps = [ps;pleGlobals.ps{j}(q,pleGlobals.q_fit)];
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
pleGlobals.mota_results = makeTable(h);

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
for j=1:length(pleGlobals.p)
    if(pleGlobals.q_fit(j))
        fprintf('p%i : %s\n', count, pleGlobals.p_labels{j});
        count = count + 1;
    end
end

%% Plot

% set relation matrix
pleGlobals.IDrelations = nan(length(pleGlobals.p));
pleGlobals.IDrelations(pleGlobals.q_fit,pleGlobals.q_fit) = ...
    pleGlobals.mota_results(1:sum(pleGlobals.q_fit),1:sum(pleGlobals.q_fit));

% cv influence
pleGlobals.IDcv = nan(1,length(pleGlobals.p));
pleGlobals.IDcv(pleGlobals.q_fit) = pleGlobals.mota_results(:,sum(pleGlobals.q_fit)+2);
% pleGlobals.IDcv(pleGlobals.q_fit) = abs((quantile(ps,0.75) - quantile(ps,0.15))./quantile(ps,0.5));

% r^2
pleGlobals.IDr2 = nan(1,length(pleGlobals.p));
pleGlobals.IDr2(pleGlobals.q_fit) = pleGlobals.mota_results(:,sum(pleGlobals.q_fit)+1);

IDweigthing = pleGlobals.IDr2; 
% IDweigthing = pleGlobals.IDcv;

tmprel = zeros([size(pleGlobals.IDrelations) 3]);
for j=1:length(pleGlobals.IDstatus)
    if(pleGlobals.q_fit(j))
        if(pleGlobals.IDstatus(j) == 1)
            tmprel(j,:,2) = pleGlobals.IDrelations(j,:).*IDweigthing(j);
        elseif(pleGlobals.IDstatus(j) == 2)
            tmprel(j,:,2) = pleGlobals.IDrelations(j,:).*IDweigthing(j);
            tmprel(j,:,1) = pleGlobals.IDrelations(j,:).*IDweigthing(j);
        elseif(pleGlobals.IDstatus(j) == 3)
            tmprel(j,:,1) = pleGlobals.IDrelations(j,:).*IDweigthing(j);
        elseif(pleGlobals.IDstatus(j) == 4)
            tmprel(j,:,3) = pleGlobals.IDrelations(j,:).*IDweigthing(j);
        end
        if(~isnan(pleGlobals.IDstatus(j)))
            tmprel(j,j,tmprel(j,j,:)~=0) = 1;
        else
            tmprel(j,j,:) = 1;
        end
    end
end
tmprel(tmprel>1) = 1;
tmprel(~pleGlobals.q_fit, :, :) = 0.5;
tmprel(:, ~pleGlobals.q_fit, :) = 0.5;

figure(figure(length(pleGlobals.p)+1))
image(tmprel);
set(gca, 'YTick', 1:length(pleGlobals.p))
set(gca, 'YTickLabel', pleGlobals.p_labels)
set(gca, 'XTick', 1:length(pleGlobals.p))
set(gca, 'XTickLabel', pleGlobals.p_labels)
%rotateticklabel(gca,45)
title('MOTA parameter relations')

% save
saveas(gcf, [pleGlobals.savePath '/mota_relations'], 'fig')
saveas(gcf, [pleGlobals.savePath '/mota_relations'], 'png')

%% save
pleSave(pleGlobals)

