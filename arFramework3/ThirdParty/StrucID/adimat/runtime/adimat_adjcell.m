% function adj = adimat_adjcell(adj, i1, i2)
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
function adj = adimat_adjcell(adj, i1, i2)
  adj = adj{i1, i2};
