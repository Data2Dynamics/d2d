% perform pleRelationsLin on Profile Likelihood Exploit
% to reveal linear related parameters
%
% pleRelationsLin

function pleRelationsLin

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

psstd = (pleGlobals.conf_ub_point(2,:) + pleGlobals.conf_lb_point(2,:))/2;

r = eye(size(ps,2));
for j=1:size(ps,2)-1
    for i=j+1:size(ps,2)
        q1 = ~isnan(pleGlobals.chi2s{i});
        q1 = q1 & transpose(pleGlobals.ps{i}(:,i)~=pleGlobals.p(i));
        a1 = mean((pleGlobals.ps{i}(q1,j) - pleGlobals.p(j))./(pleGlobals.ps{i}(q1,i) - pleGlobals.p(i)));
        
        q2 = ~isnan(pleGlobals.chi2s{j});
        q2 = q2 & transpose(pleGlobals.ps{j}(:,j)~=pleGlobals.p(j));
        a2 = mean((pleGlobals.ps{j}(q2,i) - pleGlobals.p(i))./(pleGlobals.ps{j}(q2,j) - pleGlobals.p(j)));
        
        if(psstd(j)<Inf && psstd(i)<Inf)
            a1 = a1 * (psstd(i)/psstd(j));
            a2 = a2 * (psstd(j)/psstd(i));
        elseif(psstd(j)==Inf && psstd(i)==Inf)
            a1 = 1;
            a2 = 1;
        end
        r(i,j) = (atan(a1)+atan(a2))/pi*2;
        r(i,j) = r(i,j) + 0.2*r(i,j)^3 - 0.2*r(i,j);
        
        r(j,i) = r(i,j);
    end
end

% t = r*sqrt(n-2)./sqrt(1-r.^2);
% p = tpdf(t,n);
% alpha_level = 0.01;
% qp = p<alpha_level;

%% Plot

% set relation matrix
pleGlobals.IDrelations = r.^2;

tmprel = zeros([size(pleGlobals.IDrelations) 3]);
for j=1:length(pleGlobals.IDstatus)
    if(pleGlobals.q_fit(j))
        if(pleGlobals.IDstatus(j) == 1)
            tmprel(j,:,2) = pleGlobals.IDrelations(j,:);
        elseif(pleGlobals.IDstatus(j) == 2)
            tmprel(j,:,2) = pleGlobals.IDrelations(j,:);
            tmprel(j,:,1) = pleGlobals.IDrelations(j,:);
        elseif(pleGlobals.IDstatus(j) == 3)
            tmprel(j,:,1) = pleGlobals.IDrelations(j,:);
        elseif(pleGlobals.IDstatus(j) == 4)
            tmprel(j,:,3) = pleGlobals.IDrelations(j,:);
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

% sort
dists = pdist([pleGlobals.IDrelations pleGlobals.IDrelations']);
links = linkage(dists, 'ward');
clusters = cluster(links, 'cutoff', 1);
[csort, csortindex] = sort(clusters);

figure(figure(length(pleGlobals.p)+1))
image(tmprel(csortindex,csortindex,:));
set(gca, 'YTick', 1:length(pleGlobals.p))
set(gca, 'YTickLabel', pleGlobals.p_labels(csortindex))
set(gca, 'XTick', 1:length(pleGlobals.p))
set(gca, 'XTickLabel', pleGlobals.p_labels(csortindex))
%rotateticklabel(gca,45)
title('identifiability and linear parameter relations')

% save
saveas(gcf, [pleGlobals.savePath '/lin_relations'], 'fig')
saveas(gcf, [pleGlobals.savePath '/lin_relations'], 'png')

%% save
pleSave(pleGlobals);

