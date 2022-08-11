% function [V D U]= adimat_lreigvs(A, s)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (C) 2018 Johannes Willkomm
% Copyright (C) 2015 Johannes Willkomm
function [V D U]= adimat_lreigvs(A, s)
  ws = warning('query');
  warning('off', 'MATLAB:nearlySingularMatrix');
  D = diag(s);
  k = length(s);
  U = zeros(size(A,1), k) .* A(1);
  V = zeros(k, size(A,1)) .* A(1);
  At = A.';
  for i=1:k
    ui = adimat_reigvectp(A, s(i));
    vi = adimat_reigvectp(At, s(i)).';
    vi = vi ./ (vi * ui);
    U(:,i) = ui;
    V(i,:) = vi;
  end
  warning(ws);
