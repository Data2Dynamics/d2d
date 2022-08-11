% function a_x = a_angle(a_z, x)
%
% Compute adjoint of z = angle(x), where a_z is the adjoint of z.
%
% see also a_atan2
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright (c) 2016 Johannes Willkomm
function a_x = a_angle(a_z, x)
  t1 = real(x);
  t2 = imag(x);
  [a1 a2] = a_atan2(real(a_z), 0, t2, t1);
  a_x = complex(a2, -a1);
