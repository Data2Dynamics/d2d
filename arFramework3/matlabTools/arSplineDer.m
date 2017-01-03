function ppdt = arSplineDer(pp)

% extract details from piece-wise polynomial by breaking it apart
[breaks,coefs,l,k,d] = unmkpp(pp);

% make the polynomial that describes the derivative
ppdt = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
