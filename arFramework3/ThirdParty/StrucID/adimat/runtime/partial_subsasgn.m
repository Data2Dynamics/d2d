% function J = partial_subsasgn(array, val, indices)
%
% Compute the local Jacobian of the subscripted assignment
% array(indices{:}) = val. The Jacobian is sparse logical, it will
% have exactly one 1 per column, pointing to the location of the
% component in array.
%
% see also adimat_subsasgn
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (c) 2018 Johannes Willkomm
function J = partial_subsasgn(adj, val, indices)
  testobj = zeros(size(adj));
  testobj(indices{:}) = reshape(1:numel(val), size(val));
  J = sparse(find(testobj), testobj(testobj ~= 0), 1, numel(testobj), numel(val));
