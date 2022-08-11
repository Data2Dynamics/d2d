% function str = admCamelToDash(s)
%
% Construct a string str from input s by replacing eaach upper case
% letter K with a dash and the same lower case letter -k.
%
% see also admBuildFlags
%
% This file is part of the ADiMat runtime environment
%
% Copyright 2018 Johannes Willkomm
% Copyright 2009,2011 Johannes Willkomm, Institute for Scientific Computing
%           RWTH Aachen University.
% This code is under development! Use at your own risk! Duplication,
% modification and distribution FORBIDDEN!
function str = admCamelToDash(s)
  str = '';
  upInd = find(s <= 'Z' & s >= 'A');
  offset = 1;
  for i=1:length(upInd)
    low = char('a' + s(upInd(i)) - 'A');
    str = [str s(offset:upInd(i)-1) '-' low];
    offset = upInd(i) +1;
  end
  str = [str s(offset:end)];
