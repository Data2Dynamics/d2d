% arUpdateCheckstr;
% 
% ar = arUpdateCheckstr(ar)
% 
% ar = arUpdateCheckstr(ar,dosave)
% 
%   dosave      Default: false
%               If true, then the fields used to create the checksum are
%               saved in folder Checksums
% 
%   The checksums of data-setting, parameter-setting, setup-setting and
%   fit-setting are updated.
% 
%   checkstrs.total is a cumulative checksum which is equal if all other
%   checksums are equal.

function arStruct = arUpdateCheckstr(arStruct,dosave)
global ar
if ~exist('arStruct','var') || isempty(arStruct)
    arStruct = ar;
    useglobal = true;
else
    useglobal = false;
end
if ~exist('dosave','var') || isempty(dosave)
    dosave = false;
end

% update checkstrs:
arStruct.checkstrs = struct;
arStruct.checkstrs.data = arChecksumData(arStruct,dosave);
arStruct.checkstrs.para = arChecksumPara(arStruct,dosave);
arStruct.checkstrs.fitting = arChecksumFitting(arStruct,dosave);
if isfield(arStruct,'checkstr')
    arStruct.checkstrs.fkt = arStruct.checkstr;
else
    arStruct.checkstrs.fkt = '';
end

checksum = arAddToCheckSum(arStruct.checkstrs.data);
checksum = arAddToCheckSum(arStruct.checkstrs.fitting,checksum);
checksum = arAddToCheckSum(arStruct.checkstrs.para,checksum);
arStruct.checkstrs.total = arAddToCheckSum(arStruct.checkstrs.fkt,checksum,true);
clear checksum

arStruct.checkstrs.version = '20180622'; % This number should be changed if the way checksums are calculated change and old/new ones are not comparable

if useglobal
    ar = arStruct;
end

