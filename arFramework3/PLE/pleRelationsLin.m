% perform pleRelationsLin on Profile Likelihood Exploit
% to reveal linear related parameters
%
% pleRelationsLin

function pleRelationsLin

global ar

if(~isfield(ar,'ple') || isempty(ar.ple))
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

psstd = (ar.ple.conf_ub_point(2,:) + ar.ple.conf_lb_point(2,:))/2;

r = eye(size(ps,2));
for j=1:size(ps,2)-1
    for i=j+1:size(ps,2)
        q1 = ~isnan(ar.ple.chi2s{i});
        q1 = q1 & transpose(ar.ple.ps{i}(:,i)~=ar.ple.p(i));
        a1 = mean((ar.ple.ps{i}(q1,j) - ar.ple.p(j))./(ar.ple.ps{i}(q1,i) - ar.ple.p(i)));
        
        q2 = ~isnan(ar.ple.chi2s{j});
        q2 = q2 & transpose(ar.ple.ps{j}(:,j)~=ar.ple.p(j));
        a2 = mean((ar.ple.ps{j}(q2,i) - ar.ple.p(i))./(ar.ple.ps{j}(q2,j) - ar.ple.p(j)));
        
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
ar.ple.IDrelations = r.^2;

tmprel = zeros([size(ar.ple.IDrelations) 3]);
for j=1:length(ar.ple.IDstatus)
    if(ar.qFit(j))
        if(ar.ple.IDstatus(j) == 1)
            tmprel(j,:,2) = ar.ple.IDrelations(j,:);
        elseif(ar.ple.IDstatus(j) == 2)
            tmprel(j,:,2) = ar.ple.IDrelations(j,:);
            tmprel(j,:,1) = ar.ple.IDrelations(j,:);
        elseif(ar.ple.IDstatus(j) == 3)
            tmprel(j,:,1) = ar.ple.IDrelations(j,:);
        elseif(ar.ple.IDstatus(j) == 4)
            tmprel(j,:,3) = ar.ple.IDrelations(j,:);
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

% sort
dists = pdist([ar.ple.IDrelations ar.ple.IDrelations']);
links = linkage(dists, 'ward');
clusters = cluster(links, 'cutoff', 1);
[csort, csortindex] = sort(clusters);

figure(figure(length(ar.ple.p)+1))
image(tmprel(csortindex,csortindex,:));
set(gca, 'YTick', 1:length(ar.ple.p))
set(gca, 'YTickLabel', ar.ple.p_labels(csortindex))
set(gca, 'XTick', 1:length(ar.ple.p))
set(gca, 'XTickLabel', ar.ple.p_labels(csortindex))
%rotateticklabel(gca,45)
title('identifiability and linear parameter relations')

% save
saveas(gcf, [ar.ple.savePath '/lin_relations'], 'fig')
saveas(gcf, [ar.ple.savePath '/lin_relations'], 'png')

%% save
pleSave(ar.ple);

