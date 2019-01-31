% py_struct = arToPython(struct)
%
%   Converts the ar struct to a python friendly structure.
%   
%   Essentially converts mxn cells to a 1xm cell containing m 1xn cells. 
%   There might still be data types unsupported for python import, 
%   but afaik none of them are used in the ar struct as of now.
%
%   py_struct   new python compatible ar struct
%   struct      optional struct (default = global ar)


function py_struct = arToPython(struct)

global ar

% Fill in unset optional values.
switch nargin
    case 0
        struct = ar;       
end

fields = fieldnames(struct);

py_struct = cell(1,1);

for len=1:length(struct) 
    for i=1:length(fields)
        if isstruct(struct(len).(fields{i}))
            py_struct{len}.(fields{i}) = ...
                arToPython(struct(len).(fields{i}));
        elseif iscell(struct(len).(fields{i}))
            c = size(struct(len).(fields{i}));
            if c(1) == 0
                py_struct{len}.(fields{i}) = struct(len).(fields{i});
            end
            for x=1:c(1)
                py_struct{len}.(fields{i}){x} = ...
                    struct(len).(fields{i})(x,:);
            end
        else
            py_struct{len}.(fields{i}) = struct(len).(fields{i});
        end
    end
end