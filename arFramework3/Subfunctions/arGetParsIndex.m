% p_ind = arGetParsIndex(pLabel)
%
% Get parameter index by matching label to ar.pLabel
% 
% pLabel	name of the parameter
%
% p_ind     parameter index
%
% See also arGetPars arGetParsPattern arPrint

function p_ind = arGetParsIndex(pLabel)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~iscell(pLabel))
    pLabel = {pLabel};
end

for j=1:length(pLabel)
    q = ismember(ar.pLabel, pLabel(j)); %R2013a compatible

    if(sum(q)==0)
        error('parameter %s not found', pLabel{j});
    end
    
    p_ind = find(q);

end
