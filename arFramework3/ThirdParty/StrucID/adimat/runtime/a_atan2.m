% function [a_out] = a_atan2(a_z, k, y, x)
%
% Compute adjoint of z = atan2(y, x), where a_z is the adjoint of z,
%
%   - w.r.t. y when k == 1
%   - w.r.t. x when k == 2
%
% see also a_angle
%
% Copyright (c) 2018 Johannes Willkomm
function [a_out a_out2] = a_atan2(a_z, which, y, x)
  divi = (x.^2 + y.^2);
  if any(divi(:) == 0)
    warning('atan2:both:zero', 'a_atan2(x=0,y=0)');
    divi(divi(:) == 0) = 1;
  end
  if which == 0 || which == 1
    a_y = adimat_adjred(y, a_z .* (x ./ divi));
  end
  if which == 0 || which == 2
    a_x = adimat_adjred(x, a_z .* (-y ./ divi));
  end
  switch which
   case 0
    a_out = a_y;
    a_out2 = a_x;
   case 1
    a_out = a_y;
   case 2
    a_out = a_x;
  end
