% checkstr = arChecksumPara([arStruct],[saveEvaluatedFields])
% 
% This function runs over all parameters and parameter properties like
% bounds, ar.qFit, ar.qLog10, ar.mean, ar.std and calculates a checksum.
% The parameters ar.p are NOT used for the checksum since ar.p changes
% frequently (e.g. during fitting) and the initial guess should be treated
% separately.
% 
%   arStruct        if instead of the global ar, the checksum should be
%                   evaluated for another struct, then it is provided as
%                   first argument 
%   saveEvaluatedFields  [false]
%                   if true, then a workspace is saved in folder Checksums
%                   containing the field which are evaluated for
%                   calculationg the checksum
%      

function checkstr = arChecksumPara(arStruct,saveEvaluatedFields)
global ar
if ~exist('arStruct','var') || isempty(arStruct)
    arStruct = ar;
end
if ~exist('saveEvaluatedFields','var') || isempty(saveEvaluatedFields)
    saveEvaluatedFields = false;
end



checkfields = {'qLog10','qFit','lb','ub','type','mean','std'}; % p is stored separately

if saveEvaluatedFields
    arCopy = struct;
end

checksum = [];
for i=1:length(checkfields)
    if isfield(arStruct,checkfields{i})
        val = arStruct.(checkfields{i});
        checksum = arAddToCheckSum(val,checksum);
        if saveEvaluatedFields
            arCopy.(checkfields{i}) = val;
        end
    end
end

h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

if saveEvaluatedFields
    arSaveChecksumCopy(arCopy,'para',checkstr);
end

clear checksum


