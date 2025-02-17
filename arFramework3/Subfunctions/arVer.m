function matVer = arVer

% matVer = arVer
%
% *** AVOID USING THIS FUNCTION! verLessThan IS RECOMMENDED ***
% *** DO NOT USE ver('matlab') AS THIS WON'T WORK WITH VERSIONS >= 9.10 ***
%
% Return MATLAB version with fix for 2021a (9.10) which leads to unexpected
% behaviour when converting '9.10' to double.
global ar
if isfield(ar.info,'ver')
    matVer = ar.info.ver;
else
    matVer = ver('MATLAB');
    ar.info.ver = matVer;
end

versplit = strsplit(matVer.Version,'.');
if strcmp(matVer.Version, '9.10')
    matVer.Version = '9.91';
elseif strcmp(matVer.Version, '9.11')
    matVer.Version = '9.92';
end