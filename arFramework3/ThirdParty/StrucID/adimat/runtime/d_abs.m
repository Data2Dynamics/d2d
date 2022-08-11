function [d_r r]= d_abs(d_p, p)
% D_ABS -- Compute the derivative d_p according to the value of p
%
% Copyright 2012,2013 Johannes Willkomm, Institute for Scientific Computing   
%                     TU Darmstadt
% Copyright 2001-2004 Andre Vehreschild, Institute for Scientific Computing   
%                     RWTH Aachen University
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!

  r = abs(p);

  if any(r == 0)
    warning('adimat:abs:argZero', '%s', 'd_abs(d_p, p) not defined for p==0.0');
  end

  [ndd nelx] = size(d_p);
  xr = repmat(reshape(p, [1 nelx]), [ndd 1]);
  resr = repmat(reshape(r, [1 nelx]), [ndd 1]);
  d_r = (real(xr).*real(d_p(:,:)) + imag(xr).*imag(d_p(:,:))) ./ resr;
  d_r(resr == 0) = 0;
  d_r = reshape(d_r, size(d_p));
