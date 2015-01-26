function [tUnits, response_parameter, yLabel, yNames, yUnits, iy, ...
    hys, hystds, hysss] = arGetInfo(jm, jd, jtype, linehandle_name)
global ar
tUnits = ar.model(jm).data(jd).tUnits;
response_parameter = ar.model(jm).data(jd).response_parameter;

if(jtype==1)
    yLabel = ar.model(jm).data(jd).y;
    if(isfield(ar.model(jm).data(jd), 'yNames'))
        yNames = ar.model(jm).data(jd).yNames;
    else
        yNames = [];
    end
    yUnits = ar.model(jm).data(jd).yUnits;
    iy = 1:length(ar.model(jm).data(jd).y);
    
elseif(jtype==2)
    yLabel = [ar.model(jm).u ar.model(jm).x ar.model(jm).z];
    if(isfield(ar.model(jm).data(jd), 'xNames'))
        yNames = [ar.model(jm).u ar.model(jm).xNames ar.model(jm).z];
    else
        yNames = [];
    end
    yUnits = [ar.model(jm).uUnits; ar.model(jm).xUnits; ar.model(jm).zUnits];
    iy = find([ar.model(jm).qPlotU ar.model(jm).qPlotX ar.model(jm).qPlotZ]);
    
elseif(jtype==3)
    yLabel = strrep(strrep(ar.model(jm).vs, '[', '_{'), ']', '}');
    yNames = [];
    yUnits = ar.model(jm).vUnits;
    iy = find([ar.model(jm).qPlotV]);
    
end

% get handels
if(isfield(ar.model(jm).data(jd),'plot') && ...
        isfield(ar.model(jm).data(jd).plot,linehandle_name))
    hys = ar.model(jm).data(jd).plot.(linehandle_name);
else
    hys = [];
end
if(jtype == 2 && size(ar.model(jm).data(jd).plot.(linehandle_name),1) > 1)
    hys = [];
end

if(jtype == 1 && isfield(ar.model(jm).data(jd),'plot') && ...
        isfield(ar.model(jm).data(jd).plot,'ystd'))
    hystds = ar.model(jm).data(jd).plot.ystd;
else
    hystds = [];
end
if(jtype == 2)
    if(jd~=0)
        if(isfield(ar.model(jm).data(jd),'plot') && ...
                isfield(ar.model(jm).data(jd).plot,'xss'))
            hysss = ar.model(jm).data(jd).plot.xss;
        else
            hysss = [];
        end
    else
        if(isfield(ar.model(jm).condition(jc),'plot') && ...
                isfield(ar.model(jm).condition(jc).plot,'xss'))
            hysss = ar.model(jm).condition(jc).plot.xss;
        else
            hysss = [];
        end
    end
else
    hysss = [];
end
