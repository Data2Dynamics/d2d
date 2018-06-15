% arUpdateCheckstr;
% 
% ar = arUpdateCheckstr(ar)
% 
%   The checksums of data-setting, parameter-setting, setup-setting and
%   fit-setting are updated.

function arStruct = arUpdateCheckstr(arStruct)
global ar
if nargin==0
    arStruct = ar;
end

% update checkstrs:
arStruct.checkstrs = struct;
arStruct.checkstrs.data = arChecksumData(arStruct);
arStruct.checkstrs.para = arChecksumPara(arStruct);
arStruct.checkstrs.fitting = arChecksumFitting(arStruct);
if isfield(arStruct,'checkstr')
    arStruct.checkstrs.fkt = arStruct.checkstr;
else
    arStruct.checkstrs.fkt = '';
end


