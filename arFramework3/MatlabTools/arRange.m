% y = arRange(x,dim)
% Gives range of data x in dimension dim

function y = arRange(x,dim)
if(~exist('dim','var') || isempty(dim))
    y = max(x) - min(x);
else
    y = max(x,[],dim) - min(x,[],dim); 
end
