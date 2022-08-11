% function [g_l, l]= g_eig(g_A, A)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (C) 2018 Johannes Willkomm
% Copyright (C) 2015 Johannes Willkomm
function [g_l, l]= g_eig(g_A, A)
  [V D] = eig(A);
  l = diag(D);
  g_l = sum((V \ g_A) .* V.', 2);
  g_l = reshape(g_l, size(l));
