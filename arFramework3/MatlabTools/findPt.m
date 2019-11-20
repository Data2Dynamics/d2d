%
% Small helper function to identify data points in a plot.
%
% Make sure the axis you want to find points in is active, then call this
% function and then click the point you want to know the index of.
%
function findPt
    plotted = get(gca, 'Children');
    for i = 1 : numel( plotted )
        set(plotted, 'ButtonDownFcn', @(hObject, eventdata)fnc(hObject, eventdata, i));
    end
end

function fnc(hObject, eventdata, ~)
    X = get(hObject, 'XData');
    Y = get(hObject, 'YData');
    Z = get(hObject, 'ZData');
    
    if ~isempty( Z )
        N = find( X==eventdata.IntersectionPoint(1) & Y==eventdata.IntersectionPoint(2) & Z==eventdata.IntersectionPoint(3) );
        fprintf( 'Data index %d at %g, %g, %g', N, X(N), Y(N), Z(N) );
    else
        [~,N] = min( abs(X-eventdata.IntersectionPoint(1)) + abs (Y-eventdata.IntersectionPoint(2)) );
        fprintf( 'Data index %d at %g, %g\n', N, X(N), Y(N) );
    end
    
end