function [B, i_sorted] = docluster(A, cuttoff)

% a_order=clusterdata(A,cuttoff);
% [~, i_sorted] = sort(a_order);
% B = A(i_sorted,:);

tree = linkage(A,'average');
D = pdist(A);
leafOrder = optimalleaforder(tree,D);
[~, ~, i_sorted] = dendrogram(tree,0,'Reorder',leafOrder);
close
B = A(i_sorted,:);