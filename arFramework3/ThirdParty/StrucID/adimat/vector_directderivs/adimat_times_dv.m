% function d_res = adimat_times_dv(d_A, B)
%
% Multiply derivative object d_A times B.
%
% Copyright 2012 Johannes Willkomm, Institute for Scientific Computing   
%                     TU Darmstadt
function d_res = adimat_mtimes_dv(d_A, B)
  [ndd nel] = size(d_A);
  Br = reshape(B, [1, nel]);
  Br = repmat(Br, [ndd, 1]);
  d_res = d_A(:,:) .* Br;
  d_res = reshape(d_res, size(d_A));
end

% $Id: adimat_times_dv.m 3310 2012-06-19 11:52:54Z willkomm $
