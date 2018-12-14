% arNameTrafo(str)
%
% Transforms strings to avoid unintended behavior in plot labels
% due to latex interpreter
%
%    str        string, which is transformed

function str = arNameTrafo(str)
str = strrep(str, '_', '\_');
str = strrep(str, '%', '\%');