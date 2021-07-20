function out = flatten(input) %-> [cell, struct]
    %Utility function to flatten a cell array/struct with nested cell 
    %arrays/structs.
    %
    %Parameters:
    %   input [cell, struct]:
    %       Cell array or struct to be flattened.
    %
    %Results:
    %   [cell]:
    %       Flattened cell array
    
    if ~iscell(input) && ~isstruct(input)
        error('FLATTEN:WrongInputTypeError', ...
            'Input argument must be a cell array or a struct')
    end
    
    if isstruct(input)
        input = struct2cell(input);
    end
    
    out = {};
    for i = 1:numel(input)
        if iscell(input{i}) || isstruct(input)
            out = horzcat(out, flatten(input{i}));
        elseif isstruct(input{i})
            out = horzcat(out, flatten(struct2cell(input{i})));
        else
            out = horzcat(out, input(i));
        end
    end    
end