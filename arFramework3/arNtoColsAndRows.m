% [nrows, ncols] = arNtoColsAndRows(n, rowstocols)

function [nrows, ncols] = arNtoColsAndRows(n, rowstocols)

if(~exist('rowstocols', 'var'))
    rowstocols = 0.4;
end

nrows = ceil(n^rowstocols);
if(nrows~=0)
    ncols = ceil(n / nrows);
else
    ncols = 0;
end