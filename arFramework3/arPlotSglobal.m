% Plot global sensitivities of residuals
%
% arPlotSRESglobal

function arPlotSglobal(ip, ires, iconstr)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('ip','var') || isempty(ip))
    ip = find(ar.qFit==1);
end

if(~exist('ires','var') || isempty(ires))
    ires = 1:size(ar.sres,1);
end

if(~exist('iconstr','var') || isempty(iconstr))
    iconstr = 1:size(ar.sconstr,1);
end

arChi2(true);

% constants
rowstocols = 0.5; %0.7; 0.45;
overplot = 0.1;

% rows and cols
np = length(ip);
[nrows, ncols] = arNtoColsAndRows(np, rowstocols);

% data (and prior) residuals
figure(1); clf;
ccount = 1;
for jp = ip
    g = subplot(nrows,ncols,ccount);
    
    plot(g, ar.sres(ires,jp), 'b-o');
    if(isfield(ar,'sresFD'))
        hold(g, 'on');
        plot(g, ar.sresFD(ires,jp), 'r--*');
        hold(g, 'off');
    end
    
    arSpacedAxisLimits(g, overplot);
    title(g, myNameTrafo(ar.pLabel{jp}));
    if(ccount == 1 && isfield(ar,'sresFD'))
        legend(g, {'SE','FD'});
    end
    
    if(ccount == (nrows-1)*ncols + 1)
        xlabel(g, 'residual');
        ylabel(g, 'sensitivity');
    end
    
    ccount = ccount + 1;
end

if(isfield(ar,'sresFD'))
    figure(2); clf;
    
    semilogy(min(abs(ar.sres(ires,ip) - ar.sresFD(ires,ip))./abs(ar.sres(ires,ip)), ...
        abs(ar.sres(ires,ip) - ar.sresFD(ires,ip))), 'x-')
    if(length(ip)<6)
        legend(myNameTrafo(ar.pLabel(ip)))
    end
end

% constrains
if(isempty(ar.constr))
    return
end

figure(3); clf;
ccount = 1;
for jp = ip
    g = subplot(nrows,ncols,ccount);
    
    plot(g, ar.sconstr(iconstr,jp), 'b-o');
    if(isfield(ar,'sconstrFD'))
        hold(g, 'on');
        plot(g, ar.sconstrFD(iconstr,jp), 'r--*');
        hold(g, 'off');
    end
    
    arSpacedAxisLimits(g, overplot);
    title(g, myNameTrafo(ar.pLabel{jp}));
    if(ccount == 1 && isfield(ar,'sconstrFD'))
        legend(g, {'SE','FD'});
    end
    
    if(ccount == (nrows-1)*ncols + 1)
        xlabel(g, 'constraint');
        ylabel(g, 'sensitivity');
    end
    
    ccount = ccount + 1;
end

if(isfield(ar,'sconstrFD'))
    figure(4); clf;
    
    semilogy(min(abs(ar.sconstr(iconstr,ip) - ar.sconstrFD(iconstr,ip))./abs(ar.sconstr(iconstr,ip)), ...
        abs(ar.sconstr(iconstr,ip) - ar.sconstrFD(iconstr,ip))), 'x-');
    
    if(length(ip)<6)
        legend(myNameTrafo(ar.pLabel(ip)))
    end
end





% sub-function
function str = myNameTrafo(str)
str = strrep(str, '_', '\_');
