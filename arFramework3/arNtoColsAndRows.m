function [nrows, ncols] = arNtoColsAndRows(n, rowstocols)

if(~exist('rowstocols', 'var'))
    rowstocols = 0.4;
end

nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);