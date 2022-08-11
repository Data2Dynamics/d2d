% This file is part of the ADiMat runtime environment
%
% Copyright 2011-2014 Johannes Willkomm 
%
function obj = binopFlat(obj, right, handle)
  if isobject(obj)
    isc_obj = prod(obj.m_size) == 1;
    dd1 = obj.m_derivs;
    if isobject(right)
      dd2 = right.m_derivs;
      if isc_obj
        obj.m_size = right.m_size;
      end
    else
      if isscalar(right)
        dd2 = full(right);
      else
        dd2 = reshape(full(right), [1 size(right)]);
        if isc_obj
          obj.m_size = size(right);
        end
      end
    end
  else
    dd2 = right.m_derivs;
    if isscalar(obj)
      dd1 = full(obj);
    else
      dd1 = reshape(full(obj), [1 size(obj)]);
      isc_right = prod(right.m_size) == 1;
      if isc_right
        right.m_size = size(obj);
      end
    end
    obj = right;
  end
  obj.m_derivs = bsxfun(handle, dd1, dd2);
end
