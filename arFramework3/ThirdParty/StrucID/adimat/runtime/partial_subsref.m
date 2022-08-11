% function J = partial_subsref(val, indices)
%
% Compute the local Jacobian of the subscripted selection
% val(indices{:}). The Jacobian is sparse logical, it has exactly one
% 1 per row, pointing to the origin of the component in val.
%
% see also adimat_adjsubsref
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (c) 2018 Johannes Willkomm
function J = partial_subsref(val, indices)
  testobj = reshape(1:numel(val), size(val));
  selected = testobj(indices{:});
  J = sparse(1:numel(selected), selected(:).', 1, numel(selected), numel(val));
