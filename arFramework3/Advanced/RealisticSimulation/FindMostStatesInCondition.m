function id = FindMostStatesInCondition

global ar

anzy = zeros(length(ar.model.condition),1);
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        anzy(c) = anzy(c) + sum(range(ar.model(m).condition(c).xFineSimu)>0);
    end
end
[~,id] = max(anzy);
