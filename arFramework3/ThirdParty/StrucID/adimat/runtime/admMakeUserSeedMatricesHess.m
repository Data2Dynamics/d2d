% function [seedV seedW seedRev compressedSeed coloring] = admMakeUserSeedMatricesHess(seedMatrix, admOpts)
%
% see also admHessian
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (C) 2014 Johannes Willkomm
function [seedV seedW seedRev compressedSeed coloring] = admMakeUserSeedMatricesHess(seedMatrix, admOpts)
  [seedV seedW seedRev] = admUnpackUserSeedMatricesHess(seedMatrix, admOpts);
  if ~isempty(admOpts.JPattern)
    [compressedSeed coloring] = admMakeSeedMatrixFor(seedW, admOpts.x_nCompInputs, admOpts);
  else
    compressedSeed = seedW;
    coloring = [];
  end
end
