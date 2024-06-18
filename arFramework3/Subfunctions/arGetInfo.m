% [] = arGetInfo(jm, jd, jtype, linehandle_name)
%
% Get information of data (t, names, labels, values)
%
% jype = 1: get info of obersables y
% jype = 2: get info of states x
% jype = 3: get info of v

function [tUnits, response_parameter, titles, yNames, yLabel, iy, ...
    hys, hystds, hysss] = arGetInfo(jm, jd, jtype, linehandle_name)
global ar

if(isfield(ar.model(jm), 'data'))
    tUnits = ar.model(jm).data(jd).tUnits;
    response_parameter = ar.model(jm).data(jd).response_parameter;
else
    tUnits = ar.model(jm).tUnits;
    response_parameter = 'not assigned'; % TODO
end

if(jtype==1)
    titles = ar.model(jm).data(jd).y;
    if(isfield(ar.model(jm).data(jd), 'yNames'))
        yNames = ar.model(jm).data(jd).yNames;
    else
        yNames = [];
    end
    yUnits = ar.model(jm).data(jd).yUnits;
    iy = 1:length(ar.model(jm).data(jd).y);
    
elseif(jtype==2)
    titles = [ar.model(jm).u ar.model(jm).x ar.model(jm).z];
    if(isfield(ar.model(jm), 'xNames'))
        yNames = [ar.model(jm).u ar.model(jm).xNames ar.model(jm).z];
    else
        yNames = [];
    end
    yUnits = [ar.model(jm).uUnits; ar.model(jm).xUnits; ar.model(jm).zUnits];
    iy = find([ar.model(jm).qPlotU ar.model(jm).qPlotX ar.model(jm).qPlotZ]);
    
elseif(jtype==3)
    titles = strrep(strrep(ar.model(jm).vs, '[', '_{'), ']', '}');
    if(isfield(ar.model(jm), 'v'))
        yNames = ar.model(jm).v;
    else
        yNames = [];
    end
    yUnits = ar.model(jm).vUnits;
    iy = find([ar.model(jm).qPlotV]);
end

yLabel = cell( 1, numel(iy) );
for jy = 1 : numel( iy )
    yunitetmp = '';
    if(~isempty(yUnits{iy(jy),2}))
        yunitetmp = sprintf(' [%s]', yUnits{iy(jy),2});
    end
    % Some observations may be plotted in log, but those get a special
    % legend from arGetPlotYTrafo.
    yLabel{jy} = sprintf('%s%s', yUnits{iy(jy),3}, yunitetmp);
end
    
% get handels
if(jd~=0)
    if(isfield(ar.model(jm).data(jd), 'plot') && ...
            isfield(ar.model(jm).data(jd).plot, linehandle_name))
        hys = ar.model(jm).data(jd).plot.(linehandle_name);
    else
        hys = [];
    end
else
    if(isfield(ar.model(jm), 'plot') && ...
            isfield(ar.model(jm).plot, linehandle_name))
        hys = ar.model(jm).plot.(linehandle_name);
    else
        hys = [];
    end
end

if(jd~=0)
    if(isfield(ar.model(jm).data(jd), 'plot') && ...
            isfield(ar.model(jm).data(jd).plot,'ystd'))
        hystds = ar.model(jm).data(jd).plot.ystd;
    else
        hystds = [];
    end
else
    if(isfield(ar.model(jm), 'plot') && ...
            isfield(ar.model(jm).plot,'ystd'))
        hystds = ar.model(jm).plot.ystd;
    else
        hystds = [];
    end
end

if(jd~=0 && jtype == 2)
    if(isfield(ar.model(jm).data(jd),'plot') && ...
            isfield(ar.model(jm).data(jd).plot,'xss'))
        hysss = ar.model(jm).data(jd).plot.xss;
    else
        hysss = [];
    end
else
    if(isfield(ar.model(jm),'plot') && ...
            isfield(ar.model(jm).plot,'xss'))
        hysss = ar.model(jm).plot.xss;
    else
        hysss = [];
    end
end

% Make sure that the handle direction is never ambiguous
hys = hys(:);
hystds = hystds(:);
hysss = hysss(:);