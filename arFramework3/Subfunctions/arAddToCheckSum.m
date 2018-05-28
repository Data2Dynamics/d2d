% checksum = arAddToCheckSum(val, checksum)
% 
%   This function updates a checksum according to the variable val
% 
%   checksum is an instance of class
%   java.security.MessageDigest.getInstance 
%   Initialize this object by calling this function with only a single
%   argument:
%   checksum = arAddToCheckSum(val)
% 
%   This function is not used in arCompileAll to not break compatiblity
%   with old checksums. I did not find out (quickly) why this function
%   generates other checksums. 


function checksum = arAddToCheckSum(val, checksum)
algs = {'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if(nargin<2) || isempty(checksum)
    checksum = java.security.MessageDigest.getInstance(algs{2});
end

if iscell(val) && sum(cellfun(@ischar,val)~=1)==0
    val = [val{:}];  % for compatibiltiy with old checksums
    checksum = arAddToCheckSum_core(val, checksum);
elseif iscell(val) 
    for i=1:length(val)
        checksum = arAddToCheckSum(val{i},checksum);
    end
elseif isstruct(val)
    fn = fieldnames(val);
    for i=1:length(fn)
        checksum = arAddToCheckSum(val.(fn{i}),checksum);
    end
elseif isnumeric(val)
    checksum = arAddToCheckSum_core(sprintf('%d',val),checksum);
elseif ischar(val)
    checksum = arAddToCheckSum_core(val,checksum);
elseif islogical(val)
    checksum = arAddToCheckSum_core(num2str(val),checksum);
else
    class(val)
    error('arAddToCheckSum.m: Type not yet implemented. Please do it.');
end


function checksum = arAddToCheckSum_core(str, checksum)
if(~isempty(str))
    checksum.update(uint8(str(:)));
end
