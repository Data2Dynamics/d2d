% This file is part of the ADiMat runtime environment
%
% Copyright (C) 2015,2018 Johannes Willkomm
% Copyright (C) 2011,2012,2013 Johannes Willkomm
function obj = reshape(obj, varargin)
  epos = cellfun('isempty', varargin);
  if any(epos)
    eppos = find(epos);
    varargin{eppos} = prod(obj.m_size) ./ prod(cat(1, varargin{~epos}));
  end
  obj.m_size = adimat_normalize_size([varargin{:}]);
  obj.m_derivs = reshape(full(obj.m_derivs), [obj.m_ndd obj.m_size]);
end
