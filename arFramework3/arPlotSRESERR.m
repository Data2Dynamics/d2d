% Plot models sensitivities of residuals of error model
%
% arPlotSY

function arPlotSRESERR

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

arSimu(true, true);
arSimu(true, false);

fcount = 1;
for jm = 1:length(ar.model)
    nd = length(ar.model(jm).data);
    for jd = 1:nd
        arRaiseFigure(ar.model(jm).data(jd), 'fighandel_sreserr', ['SRESERR: ' ar.model(jm).data(jd).name ' - ' ar.model(jm).data(jd).checkstr], fcount);
        
        % rows and cols
        [ncols, nrows, ny] = myColsAndRows(jm, jd);
        
        np = length(ar.model(jm).data(jd).p);
        for jy = 1:ny
            g = subplot(nrows,ncols,jy);
            arSubplotStyle(g);
            
            legendhandle = zeros(1,np);
            
            for jp = 1:np
                linestyle = arLineMarkersAndColors(jp, np, [], 'none');
                ltmp = plot(g, ar.model(jm).data(jd).tExp, ar.model(jm).data(jd).sreserr(:,jy,jp), linestyle{:}, 'Marker', 'o');
                legendhandle(jp) = ltmp;
                hold(g, 'on');
                if(isfield(ar.model(jm).data(jd), 'sresFD'))
                    plot(g, ar.model(jm).data(jd).tExp, ar.model(jm).data(jd).sreserrFD(:,jy,jp), linestyle{:}, 'Marker', '*');
                end
            end
            hold(g, 'off');
            
            arSpacedAxisLimits(g);
            title(g, arNameTrafo(ar.model(jm).data(jd).y{jy}));
            if(jy == 1)
                legend(g, legendhandle, arNameTrafo(ar.model(jm).data(jd).p));
            end
            
            if(jy == (nrows-1)*ncols + 1)
                xlabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).tUnits{3}, ar.model(jm).data(jd).tUnits{2}));
                ylabel(g, 'sensitivity');
            end
        end        
        fcount = fcount + 1;
    end
end




%% sub-functions



function [ncols, nrows, ny] = myColsAndRows(jm, jd)
global ar
ny = size(ar.model(jm).data(jd).y, 2);
[nrows, ncols] = arNtoColsAndRows(ny);


