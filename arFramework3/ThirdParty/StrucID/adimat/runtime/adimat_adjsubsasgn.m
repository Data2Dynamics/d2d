% function adj = adimat_adjsubsasgn(adj, val, varargin)
%
% Compute adjoint of subscripted assignment selection
% val(varargin{:}), where adj is the adjoint of val.
%
% see also partial_subsasgn, adimat_adjsubsref, adimat_adjred
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (c) 2018 Johannes Willkomm
function adj = adimat_adjsubsasgn(adj, val, varargin)
  if iscell(adj) || isstruct(adj)
    adj = adj(varargin{:});
  else
    J = partial_subsasgn(adj, val, varargin);
    adj = adj(:).' * J;
    adj = reshape(adj, size(val));
  end
