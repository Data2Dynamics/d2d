% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
%
function g = real(g)
  g.deriv = cellfun(@real, g.deriv, 'uniformoutput', false);
