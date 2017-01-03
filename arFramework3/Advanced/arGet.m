% data = arGet(fields, m, dset, sig, struct)
%
%   Returns data from ar struct according to a cell of field names
% 
%   data    cell array containing all fields specified in fields
%   fields  single string or cell array of strings
%   m       optional model count (default = 1)
%   dset    optional dataset(!) count (default = 1)
%   sig     optional chop ot to sig significant digits (default = false)
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
            c = 1;
            try
                data{i} = cell(1, length(struct.model(m).plot(dset).dLink));
                for j=struct.model(m).plot(dset).dLink
                    try
                        data{i}{c}=subsref(struct.model(m).data(j), substruct(S{:}));
                    catch
                        data{i}{c} = {};
                    end

                    if isnumeric(data{i}{c})
                        if length(data{i}{c}) == 1
                            data{i}{c} = {{data{i}{c}}};
                        elseif sig
                            data{i}{c} = chop(data{i}{c}, sig);
                        end
                    end
                    c = c + 1;
                end
            catch
                data{i} = {};
            end
        elseif(regexp(char(fields(i)), '^condition.'))
        
            fields(i) = {regexprep(char(fields(i)),'^condition.', '')};
            S = substruct_arg(fields(i));

            try
                data{i} = cell(1, length(struct.model(m).plot(dset).dLink));
                
                c = 1;
                for h=struct.model(m).plot(dset).dLink
                    if h == 0
                        try
                            data{i}{c}=subsref(struct.model(m).condition(1), substruct(S{:}));
                        catch
                            data{i}{c}={};
                        end
                        if isnumeric(data{i}{c})
                            if length(data{i}{c}) == 1
                                data{i}{c} = {{data{i}{c}}};
                            elseif sig
                                data{i}{c} = chop(data{i}{c}, sig);
                            end
                        end
                    else
                        for j=struct.model(m).data(h).cLink
                            try
                                data{i}{c}=subsref(struct.model(m).condition(j), substruct(S{:}));
                            catch
                                data{i}{c}={};
                            end
                            if isnumeric(data{i}{c})
                                if length(data{i}{c}) == 1
                                    data{i}{c} = {{data{i}{c}}};
                                elseif sig
                                    data{i}{c} = chop(data{i}{c}, sig);
                                end
                            end
                        c = c + 1;
                        end
                    end
                end
            catch
                data{i} = {};
            end
        elseif(regexp(char(fields(i)), '^plot.'))   
            fields(i) = {regexprep(char(fields(i)),'^plot.', '')};
            S = substruct_arg(fields(i));
            
            data{i} = cell(1, length(struct.model(m).plot));
            for j=1:length(struct.model(m).plot)
                try
                    data{i}{j}=subsref(struct.model(m).plot(j), substruct(S{:}));
                catch
                    data{i}{j}={};
                end
                if isnumeric(data{i}{j})
                    if length(data{i}{j}) == 1
                        data{i}{j} = {{data{i}{j}}};
                    end
                end
            end
        else
            S = substruct_arg(fields(i));
            try
                data{i}=subsref(struct.model(m), substruct(S{:}));
            catch
                data{i}={};
            end
            if isnumeric(data{i})
                if length(data{i}) == 1
                    data{i} = {{data{i}}};
                end
            end
        end
    else
        S = substruct_arg(fields(i));
        try
            data{i}=subsref(struct, substruct(S{:}));
        catch
            data{i}={};
        end
        if isnumeric(data{i})
            if length(data{i}) == 1
                data{i} = {{data{i}}};
            end
        end
    end
    
end
end
    
function S = substruct_arg(field)
    S = regexp(char(field), '\.', 'split');
    S(2:2:2*end)=S;
    S(1:2:end)={'.'};
end