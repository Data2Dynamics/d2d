% This function fetches the transformation function for an observable prior to plotting.
%
%   Usage:
%     trafos = arGetYTrafo(jm, jplot)
%     
%     jm    - Model index
%     jd    - Data index
%     jplot - Plot index
%       
%   Returns:
%     List of anonymous functions containing the respective transformations
%
%   Uses:
%     ar.model(jm).plot(jplot).ytrafo       - custom trafo
%     ar.model(jm).data(#).logfitting(jy)   - log10 fitting of specific observables
%     ar.model(jm).data(3).logplotting(jy)  - log10 plotting of specific observables

function trafos = arGetPlotYTrafo(jm, jd, jplot)
    global ar;

    if isfield( ar.model(jm).plot(jplot), 'ytrafo' )
        ytrafo = ar.model(jm).plot(jplot).ytrafo;
    else
        % Unit transformation. This allows us to stick to only a single code path.
        ytrafo = @(x)x;
    end
    
    % Handle all the data transformation cases.
    trafos = cell(1, size(ar.model(jm).data(jd).y, 2));
    for jy = 1 : size(ar.model(jm).data(jd).y, 2)
        if(ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = @(x) log10(ytrafo(10.^x));
        elseif(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = @(x) ytrafo(10.^x);
        elseif(~ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))                                
            trafos{jy} = @(x) log10(ytrafo(x));
        elseif(~ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = ytrafo;
        end
    end
end
