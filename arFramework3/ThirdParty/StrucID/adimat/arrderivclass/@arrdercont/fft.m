% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright 2014 Johannes Willkomm 
%
function obj = fft(obj, n, k)
  if nargin < 3
    k = adimat_first_nonsingleton(obj);
    if isscalar(obj)
      k = 2;
    end
  end
  if nargin < 2
    n = obj.m_size(k);
  end
  if admIsOctave() && k+1 > length(size(obj.m_derivs))
    % note: this concerns only trailing dimensions > 2 which are
    % singleton. hence repmat to n works
    obj.m_derivs = repmat(obj.m_derivs, adimat_repind(length(size(obj.m_derivs)), k+1,n));
  else
    obj.m_derivs = fft(obj.m_derivs, n, k+1);
  end
  obj.m_size(k) = n;
end
