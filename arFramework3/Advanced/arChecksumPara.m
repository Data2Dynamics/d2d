% This function runs over all parameters and parameter properties like
% bounds, ar.qFit, ar.qLog10, ar.mean, ar.std and calculates a checksum.
% The parameters ar.p are NOT used for the checksum since ar.p changes
% frequently (e.g. during fitting) and the initial guess should be treated
% separately.

function checkstr = arChecksumPara(arStruct)
global ar
if nargin==0
    arStruct = ar;
end
    
checkfields = {'qLog10','qFit','lb','ub','type','mean','std'}; % p is stored separately

checksum = [];
for i=1:length(checkfields)
    val = arStruct.(checkfields{i});
    checksum = arAddToCheckSum(val,checksum);
end

h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

clear checksum


