% @arrdercont/cross
% function obj = cross(a, b, dim)
%
% Compute cross product.
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (c) 2016 Johannes Willkomm
function obj = cross(a, b, dim)
  sza = size(a);
  szb = size(b);
  if isobject(a)
    if nargin < 3
      dim = find(size(a) == 3);
    end
    a.m_derivs = cross(a.m_derivs, repmat(reshape(b, [1 size(a)]), [a.m_ndd ones(1, length(size(a)))]), dim+1);
    obj = a;
  else
    if nargin < 3
      dim = find(size(b) == 3);
    end
    b.m_derivs = cross(repmat(reshape(a, [1 size(b)]), [b.m_ndd ones(1, length(size(b)))]), b.m_derivs, dim+1);
    obj = b;
  end
  if ~isequal(sza, szb)
    assert(nargin < 3);
    obj = reshape(obj, [], 3);
  end
end
