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
%     Legend for the y-axes
%
%   Uses:
%     ar.model(jm).plot(jplot).ytrafo       - custom trafo
%     ar.model(jm).data(#).logfitting(jy)   - log10 fitting of specific observables
%     ar.model(jm).data(3).logplotting(jy)  - log10 plotting of specific observables

function [trafos, ylegend] = arGetPlotYTrafo(jm, jd, jplot)
    global ar;

    % Do we have an additional transformation?
    if isfield( ar.model(jm).plot(jplot), 'ytrafo' )
        ytrafo = ar.model(jm).plot(jplot).ytrafo;
        
        % Grab textual transform representation of the function
        c = char(ytrafo);
        ffunc = strrep(strrep( c(strfind(c,')')+1:end), c(strfind(c,'(')+1:strfind(c,')')-1), '%s' ), '.', '' );
        nargs = numel( strfind( ffunc, '%s' ) );
        
        trafod = 1;
    else
        % Unit transformation. This allows us to stick to only a single code path.
        ytrafo = @(x)x;
        trafod = 0;
    end
    
    % Handle all the data transformation cases.
    trafos = cell(1, size(ar.model(jm).data(jd).y, 2));
    ylegend = cell(1, size(ar.model(jm).data(jd).y, 2));
    for jy = 1 : size(ar.model(jm).data(jd).y, 2)
        yunittmp = '';
        yUnits = ar.model(jm).data(jd).yUnits;
        obsName = yUnits{jy, 3};
        unitName = yUnits{jy, 2};
        
        % Did we transform the parameters via a plot transform? Then also
        % update the y-labels!
        if ( trafod )
            args = repmat( {obsName}, 1, nargs );
            obsName = sprintf( ffunc, args{:} );
            args = repmat( {unitName}, 1, nargs );
            unitName = sprintf( ffunc, args{:} );
            unitName = '';
        end
        
        if(~isempty(unitName))
            yunittmp = sprintf(' [%s]', unitName);
        end
        
        % Here we handle the 'regular' log10 transformations
        if(ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = @(x) log10(ytrafo(10.^x));
            ylegend{jy} = sprintf('log_{10}(%s)%s', obsName, yunittmp);
        elseif(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = @(x) ytrafo(10.^x);
            ylegend{jy} = sprintf('%s%s', obsName, yunittmp);
        elseif(~ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = @(x) log10(ytrafo(x));
            ylegend{jy} = sprintf('log_{10}(%s)%s', obsName, yunittmp);
        elseif(~ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
            trafos{jy} = ytrafo;
            ylegend{jy} = sprintf('%s%s', obsName, yunittmp);
        end
    end
end
