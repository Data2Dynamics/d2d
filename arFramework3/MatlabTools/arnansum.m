function summand = arnansum(obj, varargin)
    obj(isnan(obj)) = 0;
    summand = sum( obj, varargin{:} );
end