% This function fetches the independent variable transformation function for an observable prior to plotting.
%
%   Usage:
%     trafos = arGetXTrafo(jm, jplot)
%     
%     jm    - Model index
%     jplot - Plot index
%       
%   Returns:
%     Anonymous function containing the respective transformation
%     Legend for the x-axes
%
%   Uses:
%     ar.model(jm).plot(jplot).doseresponselog10axis       - log10?
%     ar.model(jm).plot(jplot).xtrafo                      - custom xtransform

function [xtrafo, xlegend] = arGetPlotXTrafo(jm, jplot)
    global ar;

    jd = ar.model(jm).plot(jplot).dLink(1);
    
    % Not a dose response
    if ( ar.model(jm).plot(jplot).doseresponse == 0 )
        xtrafo = @(x)x;
        if jd==0 % not plotting data
            tUnits = ar.model(jm).tUnits;            
        else
            tUnits = ar.model(jm).data(jd).tUnits;
        end
        xlegend = sprintf('%s [%s]', tUnits{3}, tUnits{2});
        return;
    end
    
    % Do we have a user specified transform?
    if(isfield(ar.model(jm).plot(jplot), 'xtrafo'))
        trafod = 1;
        xtrafo = ar.model(jm).plot(jplot).xtrafo;
        
        % Grab textual transform representation of the function
        c = char(xtrafo);
        ffunc = strrep(strrep( c(strfind(c,')')+1:end), c(strfind(c,'(')+1:strfind(c,')')-1), '%s' ), '.', '' );
        nargs = numel( strfind( ffunc, '%s' ) );
    else
        trafod = 0;
        xtrafo = @(x)x;
    end
    response_parameter = ar.model(jm).data(jd).response_parameter;
    
    if ( trafod )
        args = repmat( {response_parameter}, 1, nargs );
        response_parameter = sprintf( ffunc, args{:} );
    end
    
    if(isfield(ar.model(jm).plot(jplot), 'doseresponselog10xaxis'))
        logplotting_xaxis = ar.model(jm).plot(jplot).doseresponselog10xaxis;
        if ( logplotting_xaxis )
            xtrafo = @(x)log10( xtrafo(x) );
            xlegend = sprintf('log_{10}(%s)', arNameTrafo(response_parameter));
        else
            xlegend = sprintf('%s', arNameTrafo(response_parameter));
        end
    else
        xtrafo = @(x)log10( xtrafo(x) );
        xlegend = sprintf('log_{10}(%s)', arNameTrafo(response_parameter));
    end    
end
