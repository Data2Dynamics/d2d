% function [dpartial y] = dpartial_atan(x)
%
% Compute partial derivative diagonal of y = atan(x).
%
% see also dpartial_exp.
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright 2012 Johannes Willkomm, Fachgebiet Scientific Computing
%                     TU Darmstadt
function [dpartial y] = dpartial_atan(x)
  divi = 1 + x.^2;
  if any(divi(:) == 0)
    warning('atan:at1', 'd atan(x = +/-i)');
    divi(divi(:) == 0) = 1;
  end
  dpartial = 1 ./ divi;
  if nargout > 1
    y = atan(x);
  end
