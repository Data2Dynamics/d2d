% Plot Y models sensitivities
%
% arPlotSY

function arPlotSY

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

fcount = 1;
for jm = 1:length(ar.model)
    nd = length(ar.model(jm).data);
    for jd = 1:nd
        arRaiseFigure(ar.model(jm).plots(jf), 'fighandel_sy', ...
            ['SY: ' ar.model(jm).data(jd).name ' - ' ar.model(jm).data(jd).checkstr], fcount);
        
        % rows and cols
        ny = size(ar.model(jm).data(jd).y, 2);
        [nrows, ncols] = arNtoColsAndRows(ny);
        
        np = length(ar.model(jm).data(jd).p);
        for jy = 1:ny
            g = subplot(nrows,ncols,jy);
            arSubplotStyle(g);
            
            legendhandle = zeros(1,np);
            
            for jp = 1:np
                linestyle = arLineMarkersAndColors(jp, np, [], 'none');
                ltmp = plot(g, ar.model(jm).data(jd).tFine, ar.model(jm).data(jd).syFineSimu(:,jy,jp), linestyle{:});
                legendhandle(jp) = ltmp;
                hold(g, 'on');
                if(isfield(ar.model(jm).data(jd), 'syExpSimuFD'))
                    plot(g, ar.model(jm).data(jd).tExp, ar.model(jm).data(jd).syExpSimuFD(:,jy,jp), linestyle{:}, 'Marker', '*');
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

