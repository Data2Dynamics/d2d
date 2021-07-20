function varargout = map(handle, varargin)
    %TODO       
    
    if iscell(varargin{1})
        [varargout{1:nargout}] = cellfun(handle, varargin{:});
    elseif isnumeric(varargin{1}) || isstring(varargin{1})
        [varargout{1:nargout}] = arrayfun(handle, varargin{:});
    else
        error('MAP:WrongListTypeError', ...
            '"list" input must be of type cell or string array')
    end
end