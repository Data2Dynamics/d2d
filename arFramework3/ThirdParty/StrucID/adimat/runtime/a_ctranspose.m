% function adj = a_ctranspose(adj, val)
%
% compute adjoint of val in z = ctranspose(val), where adj is the
% adjoint of z.
%
% see also a_zeros, a_sum
%
% This file is part of the ADiMat runtime environment
%
% Copyright 2013 Johannes Willkomm
%
function adj = a_ctranspose(adj, val)
  adj = adj';
