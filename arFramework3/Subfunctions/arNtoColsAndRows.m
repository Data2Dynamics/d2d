% [nrows, ncols] = arNtoColsAndRows(n, rowstocols)
%
%   usage:
%   [nrows, ncols] = arNtoColsAndRows(n);
%   subplot(nrows, ncols, j);
%


function [nrows, ncols] = arNtoColsAndRows(n, rowstocols)

if(~exist('rowstocols', 'var'))
    rowstocols = 0.4;  % 0.5 == sqrt would means rows and cols are treated equaly, smaller numbers => less rows, larger number => more rows
end

nrows = ceil(n^rowstocols);  % ncols is rounded towards infty => more rows than cols
% nrows = floor(n^rowstocols);  % ncols is rounded towards 0 => less rows than cols
if(nrows~=0)
    ncols = ceil(n / nrows);
else
    ncols = 0;
end