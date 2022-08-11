% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright 2011-2014 Johannes Willkomm 
%
function obj = ifft(obj, varargin)
  mode = 'nonsymmetric';
  if nargin > 1 && ischar(varargin{end})
    mode = varargin{end};
    varargin = varargin(1:end-1);
  end
  if length(varargin) < 2
    k = adimat_first_nonsingleton(obj);
    if isscalar(obj)
      k = 2;
    end
  else
    k = varargin{2};
  end
  if length(varargin) < 1
    n = obj.m_size(k);
  else
    n = varargin{1};
  end
  if admIsOctave()
    modeArgs = {};
  else
    modeArgs = {mode};
  end
  if admIsOctave() && k+1 > length(size(obj.m_derivs))
    % note: this concerns only trailing dimensions > 2 which are
    % singleton. hence repmat to n works
    obj.m_derivs = repmat(obj.m_derivs, adimat_repind(length(size(obj.m_derivs)), k+1,n))./n;
  else
    obj.m_derivs = ifft(obj.m_derivs, n, k+1, modeArgs{:});
  end
  obj.m_size(k) = n;
end
