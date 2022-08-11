% function [g_V, V, g_D, D]= g_eig2(g_A, A)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (C) 2018 Johannes Willkomm
% Copyright (C) 2015 Johannes Willkomm
function [g_V, V, g_D, D] = g_eig2(g_A, A)
  [V D] = eig(A);
  TMP1 = V \ (g_A * V);
  TMP2 = d_eig_F(diag(D));
  inds = sub2ind(size(D), 1:size(D,1), 1:size(D,1));
  g_D(inds) = TMP1(inds);
  g_V = V * (TMP2 .* TMP1);
