% function [a_x a_y] = a_hypot(a_z, x, y)
%
% Compute adjoints of z = hypot(x, y), where a_z is the adjoint of z.
%
% see also a_zeros, a_mean
%
% Copyright 2018 Johannes Willkomm
% Copyright 2012,2013 Johannes Willkomm, Institute for Scientific Computing   
%                     TU Darmstadt
function [a_x, a_y] = a_hypot(a_z, x, y, noWarn)
  if nargin < 4
    noWarn = false;
  end
  t = x .^2 + y .^2;
  z = sqrt(t);
  eq0 = t == 0;
  if any(eq0(:))
    if ~noWarn
      warning('adimat:hypot:argZero', '%s', 'a_hypot(a_z, x, y) not defined for x and y both 0');
    end
    z(eq0) = 1;
  end
  a_t = a_z ./ z;
  a_x = a_t .* x;
  a_y = a_t .* y;
