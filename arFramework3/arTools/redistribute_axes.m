function h = redistribute_axes(h, n, m, isort)

hh = figure;
for j = 1:(n*m)
    subplot(n,m,j);
    plot(rand(3));
end

ichildren = [];
for j=1:length(h.Children)
    if(isa(h.Children(j), 'matlab.graphics.axis.Axes') && ~isempty(h.Children(j).Children))
        ichildren(end+1) = j; %#ok<AGROW>
    end
end

if(nargin<4)
    isort = 1:length(ichildren);
end

ccount = 1;
for j=1:length(ichildren)
%     h.Children(ichildren(j)).Position = hh.Children(n*m - isort(ccount) + 1).Position;
    h.Children(length(ichildren) - isort(ichildren(j)) + 1).Position = hh.Children(n*m - ccount + 1).Position;
    ccount = ccount + 1;
end

close(hh);