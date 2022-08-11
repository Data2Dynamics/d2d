% function [R] = adimat_leigvect(A, lambda, R)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
function R = adimat_reigvectp(A, lambda)
  X = A - lambda .* eye(size(A));
  R = complex(rand(size(A,1), 1), rand(size(A,1), 1));
  R = R ./ norm(R);
  lastStep = inf;
  tol = 1e-14;
  iter = 0;
  maxIter = 2000;
  while lastStep > tol && iter < maxIter
    T = X \ R;
    T = T ./ (norm(T) .* sign(T(1)));
    lastStep = norm(T - R);
    R = T;
    iter = iter + 1;
  end
  if iter >= maxIter
    warning('max iter %d reached, rem err %g', iter, lastStep)
  end
