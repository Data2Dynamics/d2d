% function a_A = a_eigs_1(a_l, A, varargin)
%
% Compute adjoint of A in l = eig(A), given adjoint of l.
%
% see also a_eig_10_11, a_eig_01_11, a_eig_11_11
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright 2012 Johannes Willkomm, Institute for Scientific Computing
%                     TU Darmstadt
function a_A = a_eigs_1(a_l, A, varargin)
  s = eigs(A, varargin{:});
  [V D U] = adimat_lreigvs(A, s);
  a_A = V.' * call(@diag,a_l) * U.';
