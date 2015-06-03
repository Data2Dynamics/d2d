% data = arGet(fields, jm, jd, struct)
%
%   Returns data from ar struct according to a cell of field names
% 
%   data     cell array containing all fields specified in fields
%   fields   single string or cell array of strings
%   jm       optional model count (default = 1)
%   jd       optional data count  (default = 1)
%
% Examples:
%   arGet('model.data.yFineSimu', 2, 3) 
% returns a 1x1 cell with
%   x{1} = ar.model(2).data(3).yFineSimu
%
%   arGet({'p', 'model.fy'})
% returns a 1x2 cell with
%   x{1} = ar.p, x{2} = ar.model(1).fy
%
% for convenience:
%   arGet('ar.model') == arGet('model')

function data = arGet(fields, jm, jd)

global ar

% Fill in unset optional values.
switch nargin
    case 1
        jm = 1;
        jd = 1;
    case 2
        jd = 1;
end

% make single string a cell
if(~iscell(fields))
    fields = {fields};
end

% preallocate return value for better performance
data = cell(1, length(fields));
    
for i=1:length(fields)
    
    % remove 'ar.' prefix
    fields(i) = {regexprep(char(fields(i)), '^ar.', '')};
    
    if(regexp(char(fields(i)), '^model.'))
        
        fields(i) = {regexprep(char(fields(i)),'^model.', '')};
        
        if(regexp(char(fields(i)), '^data.'))
        
            fields(i) = {regexprep(char(fields(i)),'^data.', '')};
            S = substruct_arg(fields(i));
            data{i}=subsref(ar.model(jm).data(jd), substruct(S{:}));
            
        else
            S = substruct_arg(fields(i));
            data{i}=subsref(ar.model(jm), substruct(S{:}));
        end
    else
        S = substruct_arg(fields(i));
        data{i}=subsref(ar, substruct(S{:}));
    end
    
end
end
    
function S = substruct_arg(field)
    S = regexp(char(field), '\.', 'split');
    S(2:2:2*end)=S;
    S(1:2:end)={'.'};
end