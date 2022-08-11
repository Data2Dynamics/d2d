function [g_r r]= g_abs(g_p, p)
% G_ABS -- Compute the derivative g_p according to the value of p
% g_r is the negative of g_p if p is less than zero,
% g_r is g_p if p is greater than zero, and
% if p is zero then an error is issued, because abs is not differentiable
% at zero.
%
% Copyright 2018 Johannes Willkomm
% Copyright 2012,2013 Johannes Willkomm, Institute for Scientific Computing   
%                     TU Darmstadt
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

  r = abs(p);
  eq0 = r == 0;
  divi = r;
  if any(eq0(:))
    warning('adimat:abs:argZero', '%s', 'g_abs(g_p, p) not defined for p==0.0');
    divi(eq0) = 1;
  end
  g_r = (real(p) .* call(@real, g_p) + imag(p) .* call(@imag, g_p)) ./ divi;
