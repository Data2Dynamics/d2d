% arEqualAxis
%
% makes y axis equal for the same observable names across different open 
% windows using min and max over all axis

function arEqualAxis

h = get(0,'Children');

labels = {};
lims = [];
for j=1:length(h)
    for jj=1:length(h(j).Children)
        if(strcmp(h(j).Children(jj).Type, 'axes'))
            name = h(j).Children(jj).Title.String;
            q = ismember(labels, name);
            if(sum(q)==0)
                labels{end+1} = name; %#ok<AGROW>
                lims(end+1,:) = h(j).Children(jj).YLim; %#ok<AGROW>
            else
                limstmp = h(j).Children(jj).YLim;
                lims(q,1) = min(lims(q,1), limstmp(1)); %#ok<AGROW>
                lims(q,2) = max(lims(q,2), limstmp(2)); %#ok<AGROW>
            end
        end
    end
end

for j=1:length(h)
    for jj=1:length(h(j).Children)
        if(strcmp(h(j).Children(jj).Type, 'axes'))
            name = h(j).Children(jj).Title.String;
            q = ismember(labels, name);
            
            h(j).Children(jj).YLim = lims(q,:);
        end
    end
end