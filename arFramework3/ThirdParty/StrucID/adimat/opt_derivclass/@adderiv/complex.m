% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
%
function obj = complex(a, b)
  obj = adderiv(a);
  if isscalar(a)
    obj.dims = size(b);
  end
  obj.deriv = cellfun(@complex, a.deriv, b.deriv, 'uniformoutput', false);

