% function [dpartial y] = dpartial_acos(x)
%
% Compute partial derivative diagonal of y = acos(x).
%
% see also dpartial_exp.
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright 2012 Johannes Willkomm, Fachgebiet Scientific Computing
%                     TU Darmstadt
function [dpartial y] = dpartial_acos(x)
  divi = sqrt(1 - x.^2);
  if any(divi(:) == 0)
    warning('acos:at1', 'd acos(x=1)');
    divi(divi(:) == 0) = 1;
  end
  dpartial = -1 ./ divi;
  if nargout > 1
    y = acos(x);
  end
