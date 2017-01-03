function rescale_plot(h, scale, offset)

for j=1:length(h.Children)
    if(isprop(h.Children(j), 'YData'))
        h.Children(j).YData = h.Children(j).YData*scale + offset;
    end
end