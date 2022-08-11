% function [versionNumber versionString version revision] = adimat_version(which)
%   get the version of ADiMat
%
% Copyright 2018 Johannes Willkomm
% Copyright 2010 Johannes Willkomm, Institute for Scientific Computing   
%                     RWTH Aachen University
%
% This file is part of the ADiMat runtime environment
%
function [a b c d e] = adimat_version(which)
  if nargin > 0
    [res{1:5}] = adimat_version();
    a = res{which};
  else
    a = 0.66;
    b = 'ADiMat 0.6.6-5529 (d5179372)';
    c = '0.6.6';
    d = '5529 (d5179372)';
    try
      e = d(1:min(find(d == ' '))-1);
    catch
      e = 'invalid';
    end
  end

% Local Variables:
% mode: MATLAB
% End:
