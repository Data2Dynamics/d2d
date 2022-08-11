% function adj = adimat_adjsubsref(val, adj, varargin)
%
% Compute adjoint of val, where adj is the adjoint of the subscripted
% selection val(varargin{:}).
%
% see also partial_subsref, adimat_adjsubsasgn, adimat_adjred
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (c) 2018 Johannes Willkomm
function adj = adimat_adjsubsref(val, adj, varargin)
  J = partial_subsref(val, varargin);
  adj = adj(:).' * J;
  adj = reshape(adj, size(val));
