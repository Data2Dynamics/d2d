function [B, i_sorted] = docluster(A, cutoff)

% a_order=clusterdata(A,cutoff);
% [~, i_sorted] = sort(a_order);
% B = A(i_sorted,:);

% tree = linkage(A,'average');
tree = linkage(A,'ward','euclidean','savememory','on');
D = pdist(A);
leafOrder = optimalleaforder(tree,D);
figure()
[~, ~, i_sorted] = dendrogram(tree,0,'Reorder',leafOrder);
close
B = A(i_sorted,:);