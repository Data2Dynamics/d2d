% summand = arnansum(obj, varargin)
% Count nans in object

function summand = arnansum(obj, varargin)
    obj(isnan(obj)) = 0;
    summand = sum( obj, varargin{:} );
end