function h = redistribute_axes(h, n, m)

hh = figure;
for j = 1:(n*m)
    subplot(n,m,j);
    plot(rand(3));
end

ccount = 1;
for j=1:length(h.Children)
    if(isa(h.Children(j), 'matlab.graphics.axis.Axes') && ~isempty(h.Children(j).Children))
        h.Children(j).Position = hh.Children(ccount).Position;
        ccount = ccount + 1;
    end
end

close(hh);