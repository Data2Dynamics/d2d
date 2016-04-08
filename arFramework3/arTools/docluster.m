% use a cluster function, row-wise

function [B, i_sorted] = docluster(A)

% a_order=clusterdata(A,cutoff);
% [~, i_sorted] = sort(a_order);
% B = A(i_sorted,:);

% tree = linkage(A,'average');
tree = linkage(A,'ward','euclidean','savememory','on');
if(size(A,1)<500)
    D = pdist(A);
    leafOrder = optimalleaforder(tree,D);
else
    leafOrder = [];
    warning('skipping optimalleaforder');
end

figure()
if(isempty(leafOrder))
    [~, ~, i_sorted] = dendrogram(tree,0);
else
    [~, ~, i_sorted] = dendrogram(tree,0,'Reorder',leafOrder);
end
close
B = A(i_sorted,:);