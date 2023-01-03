%  Copies the first field level of ar
function out = arCopyFit(ar)

fn = fieldnames(ar);
out = struct;
for f=1:length(fn)
    if ~isstruct(ar.(fn{f}));
        out.(fn{f}) = ar.(fn{f});
    end
end