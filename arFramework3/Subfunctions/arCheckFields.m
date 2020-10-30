function out = arCheckFields(struct, fields, silent)

% out = arCheckFields(struct, fields, [silent])
%
% Checks the presence of fields in struct
%
%   struct          struct to check
%   fields          cell array of fieldnames
%   silent          change response behavior [0]
%                       0: print message for success and error for failure
%                       1: do not print anythign or throw an error
%                       2: only throw error when not all fields are present
%
%   out             ell array of missing fieldnames

if ~exist('silent') || isempty(silent)
    silent = 0;
end

out = {};
fields = cellstr(fields);
n = length(fields);
for i = 1:n
    subfields = strsplit(fields{i}, '.');
    if numel(subfields) == 1
        isf = isfield(struct, subfields);
        missing_names = subfields(~isf);
        if ~isempty(missing_names)
            out(end+1) = missing_names;
        end
    else
        out1 = arCheckFields(struct.(subfields{1}), strjoin(subfields(2:end),'.'),1);
        if ~isempty(out1)
            out(end+1) = out1;
        end
    end
end

if ~isempty(out) && silent ~= 1
    error('arCheckFields: The following fields in are missing: %s', strjoin(out, ', '))
end
if isempty(out) && silent == 0
    fprintf('arCheckFields: All fields present\n')
end
end