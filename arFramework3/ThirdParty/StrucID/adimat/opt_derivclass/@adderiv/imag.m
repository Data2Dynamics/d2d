% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
%
function g = imag(g)
  g.deriv = cellfun(@imag, g.deriv, 'uniformoutput', false);
