% data = arGet(fields, m, dset, sig, struct)
%
%   Returns data from ar struct according to a cell of field names
% 
%   data     cell array containing all fields specified in fields
%   fields   single string or cell array of strings
%   m       optional model count (default = 1)
%   d       optional data count  (default = 1)
%   c       optional condition count (default = 1)
%   struct  optional struct (with ar structure) (default = global ar)
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

function data = arGet(fields, m, dset, sig, struct)

global ar

% Fill in unset optional values.
switch nargin
    case 1
        m = 1;
        dset = 1;
        sig = false;
        struct = ar;
    case 2
        dset = 1;
        sig = false;
        struct = ar;
    case 3
        sig = false;
        struct = ar;
    case 4
        struct = ar;
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
            
            data{i} = cell(1, length(struct.model(m).plot(dset).dLink));
            for j=struct.model(m).plot(dset).dLink
                data{i}{j}=subsref(struct.model(m).data(j), substruct(S{:}));
                
                if isnumeric(data{i}{j}) && sig
                    
                    data{i}{j} = chop(data{i}{j}, sig);
                end
            end
            
        elseif(regexp(char(fields(i)), '^condition.'))
        
            fields(i) = {regexprep(char(fields(i)),'^condition.', '')};
            S = substruct_arg(fields(i));
            
            data{i} = cell(1, length(struct.model(m).condition));
            for j=1:length(struct.model(m).condition)
                data{i}{j}=subsref(struct.model(m).condition(j), substruct(S{:}));
                
                if isnumeric(data{i}{j}) && sig
                    data{i}{j} = chop(data{i}{j}, sig);
                end
            end
        elseif(regexp(char(fields(i)), '^plot.'))   
            fields(i) = {regexprep(char(fields(i)),'^plot.', '')};
            S = substruct_arg(fields(i));
            
            data{i} = cell(1, length(struct.model(m).plot));
            for j=1:length(struct.model(m).plot)
                data{i}{j}=subsref(struct.model(m).plot(j), substruct(S{:}));
            end
        else
            S = substruct_arg(fields(i));
            data{i}=subsref(struct.model(m), substruct(S{:}));
        end
    else
        S = substruct_arg(fields(i));
        data{i}=subsref(struct, substruct(S{:}));
    end
    
end
end
    
function S = substruct_arg(field)
    S = regexp(char(field), '\.', 'split');
    S(2:2:2*end)=S;
    S(1:2:end)={'.'};
end