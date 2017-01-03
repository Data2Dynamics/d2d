function arRaiseFigs

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

for m = 1:length(ar.model)
    for j = 1:length(ar.model(m).plot)
        if(isfield(ar.model(m).plot(j),'fighandel_y'))
            figure(ar.model(m).plot(j).fighandel_y);
        end
        if(isfield(ar.model(m).plot(j),'fighandel_x'))
            figure(ar.model(m).plot(j).fighandel_x);
        end
        if(isfield(ar.model(m).plot(j),'fighandel_v'))
            figure(ar.model(m).plot(j).fighandel_v);
        end
    end
end
