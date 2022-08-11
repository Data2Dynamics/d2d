% function [L] = adimat_leigvect(A, lambda, R)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
function L = adimat_leigvectp(A, lambda)
  L = adimat_reigvectp(A.', lambda).';
