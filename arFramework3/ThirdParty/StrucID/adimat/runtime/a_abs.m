% function a_x = a_abs(a_z, x)
%
% Compute adjoint of z = abs(x), where a_z is the adjoint of z.
%
% see also a_zeros, a_mean
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright (c) 2017 Johannes Willkomm
% Copyright (c) 2016 Johannes Willkomm
% Copyright (c) 2012 Johannes Willkomm, Institute for Scientific Computing   
%                     TU Darmstadt

function a_x = a_abs(a_z, x)
  if any(x == 0)
    warning('adimat:abs:argZero', '%s', 'a_abs(a_z, x) not defined for x==0.0');
  end
  [a_r a_i] = a_hypot(call(@real, a_z), real(x), imag(x), true);
  a_x = complex(a_r, -a_i);
