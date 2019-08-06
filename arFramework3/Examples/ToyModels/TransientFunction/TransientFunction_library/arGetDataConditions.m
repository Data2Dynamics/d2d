function dc = arGetDataConditions

global ar

dc = struct;
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        if ~isempty(ar.model(m).data(d).condition)
            pat = ar.model(m).data(d).condition.parameter;
            rep = ar.model(m).data(d).condition.value;
            
            if ~isfield(dc,pat)
                dc.(pat) = {rep};
            else
                dc.(pat){end+1} = rep;
            end
        end
    end
end

fn = fieldnames(dc);
for f=1:length(fn)
    dc.(fn{f}) = unique(dc.(fn{f}));
end
