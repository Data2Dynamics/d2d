% function [V D U] = adimat_lreigs(A, varargin)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (C) 2018 Johannes Willkomm
% Copyright (C) 2015 Johannes Willkomm
function [V D U]= adimat_lreigs(A, varargin)
  ws = warning('query');
  warning('off', 'MATLAB:nearlySingularMatrix');
  [U D] = eigs(A, varargin{:});
  k = size(U, 2);
  V = zeros(size(U,2), size(U,1));
  At = A.';
  for i=1:k
    V(i,:) = adimat_reigvectp(At, D(i,i)).';
    V(i,:) = V(i,:) ./ abs(V(i,:) * U(:,i));
  end
  warning(ws);
