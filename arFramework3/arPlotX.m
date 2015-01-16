
% Plot models X
%
% arPlotX(saveToFile, fastPlot)
%
% saveToFile    [false]
% fastPlot      [false]
% 
%   After clicking the subplot of interest, the following command provides
%   annotation of the displayed plot:
%   get(gca,'UserData') 

function arPlotX(saveToFile, fastPlot)

matVer = ver('MATLAB');

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('saveToFile','var'))
    saveToFile = false;
end
if(~exist('fastPlot','var'))
    fastPlot = false;
end

% constants
if(ar.config.ploterrors == -1)
    linesize = 0.5;
else
    linesize = 2;
end

if(isfield(ar.config,'nfine_dr_plot'))
    nfine_dr_plot = ar.config.nfine_dr_plot;
    nfine_dr_method = ar.config.nfine_dr_method;
else
    nfine_dr_plot = 1;
    nfine_dr_method = 'spline';
end


if(isfield(ar.config,'plot_x_collected'))
    plot_x_collected = ar.config.plot_x_collected;
else
    plot_x_collected = false;
end

clinks = cell(size(ar.model));

figcount = 1;
for jm = 1:length(ar.model)
    clinks{jm} = cell(size(ar.model(jm).plot));
    for jplot = 1:length(ar.model(jm).plot)
        if(ar.model(jm).qPlotXs(jplot)==1)
            if(ar.config.ploterrors == -1)
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['CI-X: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
            else
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['X: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
            end
            
            if(isfield(ar.model(jm).plot(jplot), 'doseresponselog10xaxis'))
                logplotting_xaxis = ar.model(jm).plot(jplot).doseresponselog10xaxis;
            else
                logplotting_xaxis = true;
            end
            
            % plotting
            ccount = 1;
            if(~ar.model(jm).plot(jplot).doseresponse)
                % rows and cols
                [ncols, nrows, nu, nx, ~, iu, ix, iz] = myColsAndRows(jm);
                
                if(plot_x_collected)
                    cclegendstyles = zeros(1, length(ar.model(jm).u) + length(ar.model(jm).x) + length(ar.model(jm).z));
                else
                    cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                end
                
                for jd = ar.model(jm).plot(jplot).dLink
                    [t, u, x, z, ulb, uub, xlb, xub, zlb, zub, jc, dxdt] = getData(jm, jd);
                    
                    clinks{jm}{jplot} = [clinks{jm}{jplot},jc];
                    
                    Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount);
                    
                    countu = 0;
                    for ju = iu
                        countu = countu + 1;
                        if(plot_x_collected)
                            Clines = arLineMarkersAndColors(countu, length(ar.model(jm).u) + length(ar.model(jm).x) + ...
                                length(ar.model(jm).z), [], 'none');
                        end
                        if(~fastPlotTmp)
                            if(plot_x_collected)
                                g = subplot(1,1,1);
                            else
                                g = subplot(nrows,ncols,countu);
                            end
                            ar.model(jm).plot(jplot).gu(ju) = g;
                            arSubplotStyle(g);
                            ltmp = plot(g, t, u(:,ju), Clines{:}, 'LineWidth', linesize);
                            if(plot_x_collected)
                                cclegendstyles(countu) = ltmp;
                            else
                                cclegendstyles(ccount) = ltmp;
                            end
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.u(ju,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.u(ju,jc) = ltmp;
                            end
                            hold(g, 'on');
                            if(ar.config.ploterrors == -1)
                                tmpx = [t(:); flipud(t(:))];
                                tmpy = [uub(:,ju); flipud(ulb(:,ju))];
                                ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                            end
                        else
                            if(jd~=0)
                                set(ar.model(jm).data(jd).plot.u(ju,jc), 'YData', u(:,ju));
                            else
                                set(ar.model(jm).condition(jc).plot.u(ju,jc), 'YData', u(:,ju));
                            end
                        end
                    end
                    countx = 0;
                    for jx = ix
                        countx = countx + 1;
                        if(plot_x_collected)
                            Clines = arLineMarkersAndColors(countu+countx, length(ar.model(jm).u) + length(ar.model(jm).x) + ...
                                length(ar.model(jm).z), [], 'none');
                        end
                        if(~fastPlotTmp)
                            if(plot_x_collected)
                                g = subplot(1,1,1);
                            else
                                g = subplot(nrows,ncols,countx+nu);
                            end
                            ar.model(jm).plot(jplot).gx(jx) = g;
                            arSubplotStyle(g);
                            
                            % plot ssa
                            if(isfield(ar.model(jm).condition(jc), 'xFineSSA'))
                                for jssa = 1:size(ar.model(jm).condition(jc).xFineSSA_lb, 3)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ar.model(jm).condition(jc).xFineSSA_ub(:,jx,jssa); ...
                                        flipud(ar.model(jm).condition(jc).xFineSSA_lb(:,jx,jssa))];
                                    patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.2+0.8, 'FaceColor', Clines{2}*0.2+0.8)
                                    %                                 patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.4+0.6, 'FaceColor', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                for jssa = 1:size(ar.model(jm).condition(jc).xFineSSA_lb, 3)
                                    plot(t, ar.model(jm).condition(jc).xFineSSA(:,jx,jssa), 'Color', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                if(size(ar.model(jm).condition(jc).xFineSSA,3)>1)
                                    plot(t, mean(ar.model(jm).condition(jc).xFineSSA(:,jx,:),3), '--', 'Color', Clines{2})
                                end
                                
                                % density plot
                                %                             xtmp = linspace(min(min(ar.model(jm).condition(jc).xFineSSA_lb(:,jx,:)))*0.8, ...
                                %                                 max(max(ar.model(jm).condition(jc).xFineSSA_ub(:,jx,:)))*1.2, 100);
                                %                             xSSAdens = zeros(length(t), length(xtmp));
                                %                             for jssa=1:length(t)
                                %                                 xSSAdens(jssa,:) = ksdensity(squeeze((ar.model(jm).condition(jc).xFineSSA_ub(jssa,jx,:) + ...
                                %                                 ar.model(jm).condition(jc).xFineSSA_lb(jssa,jx,:))/2), xtmp);
                                %                                 xSSAdens(jssa,:) = xSSAdens(jssa,:) / sum(xSSAdens(jssa,:));
                                %                             end
                                %                             [A,B] = meshgrid(xtmp,t);
                                %                             colormap jet
                                %                             contourf(B, A, xSSAdens);
                                %                             hold(g, 'on');
                            end
                            
                            ltmp = plot(g, t, x(:,jx), Clines{:}, 'LineWidth', linesize);
                            if(plot_x_collected)
                                cclegendstyles(countu+countx) = ltmp;
                            else
                                cclegendstyles(ccount) = ltmp;
                            end
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.x(jx,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.x(jx,jc) = ltmp;
                            end
                            hold(g, 'on');
                            
                            % steady state
                            xss = x(1,jx) + dxdt(jx)*(t-min(t));
                            ltmp = plot(g, t, xss, '--', Clines{:});
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.xss(jx,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.xss(jx,jc) = ltmp;
                            end
                            
                            if(ar.config.ploterrors == -1)
                                tmpx = [t(:); flipud(t(:))];
                                tmpy = [xub(:,jx); flipud(xlb(:,jx))];
                                ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                            end
                        else
                            % steady state
                            xss = x(1,jx) + dxdt(jx)*(t-min(t));
                            xss(xss<0) = nan;
                            xss(xss>2*max(x(:,jx))) = nan;
                            
                            if(jd~=0)
                                set(ar.model(jm).data(jd).plot.x(jx,jc), 'YData', x(:,jx));
                                set(ar.model(jm).data(jd).plot.xss(jx,jc), 'YData', xss);
                            else
                                set(ar.model(jm).condition(jc).plot.x(jx,jc), 'YData', x(:,jx));
                                set(ar.model(jm).condition(jc).plot.xss(jx,jc), 'YData', xss);
                            end
                        end
                    end
                    countz = 0;
                    for jz = iz
                        countz = countz + 1;
                        if(plot_x_collected)
                            Clines = arLineMarkersAndColors(countu+countx+countz, length(ar.model(jm).u) + length(ar.model(jm).x) + ...
                            length(ar.model(jm).z), [], 'none');
                        end
                        if(~fastPlotTmp)
                            if(plot_x_collected)
                                g = subplot(1,1,1);
                            else
                                g = subplot(nrows,ncols,countz+nu+nx);
                            end
                            ar.model(jm).plot(jplot).gz(jz) = g;
                            arSubplotStyle(g);
                            
                            % plot ssa
                            if(isfield(ar.model(jm).condition(jc), 'zFineSSA'))
                                for jssa = 1:size(ar.model(jm).condition(jc).zFineSSA_lb, 3)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ar.model(jm).condition(jc).zFineSSA_ub(:,jz,jssa); ...
                                        flipud(ar.model(jm).condition(jc).zFineSSA_lb(:,jz,jssa))];
                                    patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.2+0.8, 'FaceColor', Clines{2}*0.2+0.8)
                                    %                                 patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.4+0.6, 'FaceColor', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                for jssa = 1:size(ar.model(jm).condition(jc).zFineSSA_lb, 3)
                                    plot(t, ar.model(jm).condition(jc).zFineSSA(:,jz,jssa), 'Color', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                if(size(ar.model(jm).condition(jc).zFineSSA,3)>1)
                                    plot(t, mean(ar.model(jm).condition(jc).zFineSSA(:,jz,:),3), '--', 'Color', Clines{2})
                                end
                                
                                % density plot
                                %                             xtmp = linspace(min(min(ar.model(jm).condition(jc).zFineSSA_lb(:,jz,:)))*0.8, ...
                                %                                 max(max(ar.model(jm).condition(jc).zFineSSA_ub(:,jz,:)))*1.2, 100);
                                %                             zSSAdens = zeros(length(t), length(xtmp));
                                %                             for jssa=1:length(t)
                                %                                 zSSAdens(jssa,:) = ksdensity(squeeze((ar.model(jm).condition(jc).zFineSSA_ub(jssa,jz,:) + ...
                                %                                 ar.model(jm).condition(jc).zFineSSA_lb(jssa,jz,:))/2), xtmp);
                                %                                 zSSAdens(jssa,:) = zSSAdens(jssa,:) / sum(zSSAdens(jssa,:));
                                %                             end
                                %                             [A,B] = meshgrid(xtmp,t);
                                %                             colormap jet
                                %                             contourf(B, A, zSSAdens);
                                %                             hold(g, 'on');
                            end
                            
                            ltmp = plot(g, t, z(:,jz), Clines{:}, 'LineWidth', linesize);
                            if(plot_x_collected)
                                cclegendstyles(countu+countx+countz) = ltmp;
                            else
                                cclegendstyles(ccount) = ltmp;
                            end
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.z(jz,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.z(jz,jc) = ltmp;
                            end
                            hold(g, 'on');
                            
                            if(ar.config.ploterrors == -1)
                                tmpx = [t(:); flipud(t(:))];
                                tmpy = [zub(:,jz); flipud(zlb(:,jz))];
                                ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                            end
                        else
                            if(jd~=0)
                                set(ar.model(jm).data(jd).plot.z(jz,jc), 'YData', z(:,jz));
                            else
                                set(ar.model(jm).condition(jc).plot.z(jz,jc), 'YData', z(:,jz));
                            end
                        end
                    end
                    ccount = ccount + 1;
                end
            else
                times = [];
                for jd = ar.model(jm).plot(jplot).dLink
					times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                    [ncols, nrows, nu, nx, ~, iu, ix, iz] = myColsAndRows(jm);
                end
                
                if(str2double(matVer.Version)>=8.1)
                    [conditions, iconditions, jconditions] = unique(ar.model(jm).plot(jplot).condition,'legacy'); %#ok<ASGLU>
                else
                    [conditions, iconditions, jconditions] = unique(ar.model(jm).plot(jplot).condition); %#ok<ASGLU>
                end
                
                if(plot_x_collected)
                    cclegendstyles = zeros(1, length(iu) + length(ix) + length(iz));
                else
                    cclegendstyles = zeros(1,length(times)*length(conditions));
                end
                
                for jt = 1:length(times)
                    if(isempty(conditions))
                        jcs = 1;
                    else
                        jcs = 1:length(conditions);
                    end
                    for jc = jcs
                        if(isempty(conditions))
                            ds = ar.model(jm).plot(jplot).dLink;
                        else
                            ds = ar.model(jm).plot(jplot).dLink(find(jconditions==jc)); %#ok<FNDSB>
                        end
                        
                        jd = ds(1);
                        
                        Clines = myLineStyle(length(times)*length(jcs), ccount);
                        
                        countu = 0;
                        for ju = iu
                            countu = countu + 1;
                            if(plot_x_collected)
                                Clines = arLineMarkersAndColors(countu, length(ar.model(jm).u) + length(ar.model(jm).x) + ...
                                    length(ar.model(jm).z), [], 'none');
                            end
                            [t, u, lb, ub, zero_break] = getDataDoseResponseU(jm, ju, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                            if(length(unique(t))==1)
                                t = [t-0.1; t+0.1];
                                u = [u; u]; %#ok<AGROW>
                                lb = [lb; lb]; %#ok<AGROW>
                                ub = [ub; ub]; %#ok<AGROW>
                            elseif(nfine_dr_plot>10)
                                tf = linspace(min(t), max(t), nfine_dr_plot);
                                [t, qit] = unique(t);
                                u = u(qit);
                                u = interp1(t,u,tf,nfine_dr_method);
                                if(~isempty(lb))
                                    lb = lb(qit);
                                    ub = ub(qit);
                                    lb = interp1(t,lb,tf,nfine_dr_method);
                                    ub = interp1(t,ub,tf,nfine_dr_method);
                                end
                                t = tf;
                            end
                            
                            if(~fastPlotTmp)
                                if(plot_x_collected)
                                    g = subplot(1,1,1);
                                else
                                    g = subplot(nrows,ncols,countu);
                                end
                                ar.model(jm).plot(jplot).gu(ju) = g;
                                arSubplotStyle(g);
                                ltmp = plot(g, t, u, Clines{:}, 'LineWidth', linesize);
                                if(plot_x_collected)
                                    cclegendstyles(countu) = ltmp;
                                else
                                    cclegendstyles(ccount) = ltmp;
                                end
                                ar.model(jm).data(jd).plot.u(ju,jt,jc) = ltmp;
                                hold(g, 'on');
                                if(ar.config.ploterrors == -1)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ub(:); flipud(lb(:))];
                                    ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                end
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                set(ar.model(jm).data(jd).plot.u(ju,jt,jc), 'YData', u);
                            end
                        end
                        countx = 0;
                        for jx = ix
                            countx = countx + 1;
                            if(plot_x_collected)
                                Clines = arLineMarkersAndColors(countu+countx, length(ar.model(jm).u) + length(ar.model(jm).x) + ...
                                    length(ar.model(jm).z), [], 'none');
                            end
                            [t, x, lb, ub, zero_break] = getDataDoseResponseX(jm, jx, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                            if(length(unique(t))==1)
                                t = [t-0.1; t+0.1];
                                x = [x; x]; %#ok<AGROW>
                                lb = [lb; lb]; %#ok<AGROW>
                                ub = [ub; ub]; %#ok<AGROW>
                            elseif(nfine_dr_plot>10)
                                tf = linspace(min(t), max(t), nfine_dr_plot);
                                [t, qit] = unique(t);
                                x = x(qit);
                                x = interp1(t,x,tf,nfine_dr_method);
                                if(~isempty(lb))
                                    lb = lb(qit);
                                    ub = ub(qit);
                                    lb = interp1(t,lb,tf,nfine_dr_method);
                                    ub = interp1(t,ub,tf,nfine_dr_method);
                                end
                                t = tf;
                            end
                            
                            if(~fastPlotTmp)
                                if(plot_x_collected)
                                    g = subplot(1,1,1);
                                else
                                    g = subplot(nrows,ncols,countx+nu);
                                end
                                ar.model(jm).plot(jplot).gx(jx) = g;
                                arSubplotStyle(g);
                                ltmp = plot(g, t, x, Clines{:}, 'LineWidth', linesize);
                                if(plot_x_collected)
                                    cclegendstyles(countu+countx) = ltmp;
                                else
                                    cclegendstyles(ccount) = ltmp;
                                end
                                ar.model(jm).data(jd).plot.x(jx,jt,jc) = ltmp;
                                hold(g, 'on');
                                if(ar.config.ploterrors == -1)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ub(:); flipud(lb(:))];
                                    ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                end
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                set(ar.model(jm).data(jd).plot.x(jx,jt,jc), 'YData', x);
                            end
                        end
                        countz = 0;
                        for jz = iz
                            countz = countz + 1;
                            if(plot_x_collected)
                                Clines = arLineMarkersAndColors(countu+countx+countz, length(ar.model(jm).u) + length(ar.model(jm).x) + ...
                                    length(ar.model(jm).z), [], 'none');
                            end
                            [t, z, lb, ub, zero_break] = getDataDoseResponseZ(jm, jz, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                            if(length(unique(t))==1)
                                t = [t-0.1; t+0.1];
                                z = [z; z]; %#ok<AGROW>
                                lb = [lb; lb]; %#ok<AGROW>
                                ub = [ub; ub]; %#ok<AGROW>
                            elseif(nfine_dr_plot>10)
                                tf = linspace(min(t), max(t), nfine_dr_plot);
                                [t, qit] = unique(t);
                                z = z(qit);
                                z = interp1(t,z,tf,nfine_dr_method);
                                if(~isempty(lb))
                                    lb = lb(qit);
                                    ub = ub(qit);
                                    lb = interp1(t,lb,tf,nfine_dr_method);
                                    ub = interp1(t,ub,tf,nfine_dr_method);
                                end
                                t = tf;
                            end
                            
                            if(~fastPlotTmp)
                                if(plot_x_collected)
                                    g = subplot(1,1,1);
                                else
                                    g = subplot(nrows,ncols,countz+nu+nx);
                                end
                                ar.model(jm).plot(jplot).gz(jz) = g;
                                arSubplotStyle(g);
                                ltmp = plot(g, t, z, Clines{:}, 'LineWidth', linesize);
                                if(plot_x_collected)
                                    cclegendstyles(countu+countx+countz) = ltmp;
                                else
                                    cclegendstyles(ccount) = ltmp;
                                end
                                ar.model(jm).data(jd).plot.z(jz,jt,jc) = ltmp;
                                hold(g, 'on');
                                if(ar.config.ploterrors == -1)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ub(:); flipud(lb(:))];
                                    ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                end
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                set(ar.model(jm).data(jd).plot.z(jz,jt,jc), 'YData', z);
                            end
                        end
                        ccount = ccount + 1;
                    end
                end
            end
            
            % axis & titles
            if(~fastPlotTmp && exist('suptitle','file')==2) % suptitle function is available (can be downloaded from matlab fileexchange)
                suptitle(arNameTrafo([ar.model(jm).name,': ',ar.model(jm).plot(jplot).name]),'FontSize',12)
            end
            
            jd = ar.model(jm).plot(jplot).dLink(1);
            if(isfield(ar.model(jm), 'data'))
                for jc = 1:length(ar.model(jm).data(jd).condition)
                    if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
                        jcondi = jc;
                    end
                end
            end
            countu = 0;
            for ju = iu
                countu = countu + 1;
                g = ar.model(jm).plot(jplot).gu(ju);
                
                set(g,'UserData',...
                    struct('jm',jm,'jplot',jplot,'ju',ju, ...
                    'dLink',ar.model(jm).plot(jplot).dLink, ...
                    'cLink',clinks{jm}{jplot}, ...
                    'u', ar.model(jm).u{ju},...
                    'model_name',ar.model(jm).name, ...
                    'plot_name',ar.model(jm).plot(jplot).name ...
                    ))
                
                if(~fastPlotTmp)
                    hold(g, 'off');
                    
                    title(g, arNameTrafo(ar.model(jm).u{ju}));
                    if(ju == iu(1))
                        if(plot_x_collected)
                            legend(g, cclegendstyles, arNameTrafo([ar.model(jm).u(iu) ar.model(jm).x(ix) ar.model(jm).z(iz)]))
                        else
                            if((~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                                if(~ar.model(jm).plot(jplot).doseresponse)
                                    if(length(ar.model(jm).plot(jplot).dLink)>1)
                                        legend(g, cclegendstyles, arNameTrafo(ar.model(jm).plot(jplot).condition))
                                    end
                                else
                                    legendtmp = {};
                                    ccount = 1;
                                    for jt=1:length(times)
                                        if(~isempty(conditions))
                                            for jc = 1:length(conditions)
                                                legendtmp{ccount} = sprintf('t=%g%s : %s', times(jt), ar.model(jm).tUnits{2}, conditions{jc}); %#ok<AGROW>
                                                ccount = ccount + 1;
                                            end
                                        else
                                            legendtmp{ccount} = sprintf('t=%g%s', times(jt), ar.model(jm).tUnits{2}); %#ok<AGROW>
                                            ccount = ccount + 1;
                                        end
                                    end
                                    legend(g, cclegendstyles, arNameTrafo(legendtmp))
                                end
                            end
                        end
                    end
                    if(nx==0 && countu == (nrows-1)*ncols + 1)
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                        else
                            if(logplotting_xaxis)
                                xlabel(g, sprintf('log_{10}(%s)', arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter)));
                            else
                                xlabel(g, sprintf('%s', arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter)));
                            end
                        end
                    end
                    ylabel(g, sprintf('%s [%s]', ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2}));
                end
                arSpacedAxisLimits(g);
            end
            countx = 0;
            for jx = ix
                countx = countx + 1;
                g = ar.model(jm).plot(jplot).gx(jx);
                
                set(g,'UserData',...
                    struct('jm',jm,'jplot',jplot,'jx',jx, ...
                    'dLink',ar.model(jm).plot(jplot).dLink, ...
                    'cLink',clinks{jm}{jplot}, ...
                    'xName',ar.model(jm).xNames{jx}, ...
                    'model_name',ar.model(jm).name, ...
                    'plot_name',ar.model(jm).plot(jplot).name ...
                    ))
                
                if(~fastPlotTmp)
                    hold(g, 'off');    
                    if(nu == 0 && jx == ix(1))
                        if(plot_x_collected)
                            legend(g, cclegendstyles, arNameTrafo([ar.model(jm).u(iu) ar.model(jm).x(ix) ar.model(jm).z(iz)]))
                        else
                            if((~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                                if(~ar.model(jm).plot(jplot).doseresponse)
                                    if(length(ar.model(jm).plot(jplot).dLink)>1)
                                        legend(g, cclegendstyles, arNameTrafo(ar.model(jm).plot(jplot).condition))
                                    end
                                else
                                    legendtmp = {};
                                    ccount = 1;
                                    for jt=1:length(times)
                                        if(~isempty(conditions))
                                            for jc = 1:length(conditions)
                                                legendtmp{ccount} = sprintf('t=%g%s : %s', times(jt), ar.model(jm).tUnits{2}, conditions{jc}); %#ok<AGROW>
                                                ccount = ccount + 1;
                                            end
                                        else
                                            legendtmp{ccount} = sprintf('t=%g%s', times(jt), ar.model(jm).tUnits{2}); %#ok<AGROW>
                                            ccount = ccount + 1;
                                        end
                                    end
                                    legend(g, cclegendstyles, arNameTrafo(legendtmp))
                                end
                            end
                        end
                    end

                    if(isfield(ar.model(jm), 'xNames') && ~isempty(ar.model(jm).xNames{jx}) && ...
                            ~strcmp(ar.model(jm).xNames{jx},ar.model(jm).x{jx}))
                        title(g, [arNameTrafo(ar.model(jm).xNames{jx}) ' (' arNameTrafo(ar.model(jm).x{jx}) ')']);
                    else
                        title(g, arNameTrafo(ar.model(jm).x{jx}));
                    end
                    if(countx+nu == (nrows-1)*ncols + 1)
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                        else
                            if(logplotting_xaxis)
                                xlabel(g, sprintf('log_{10}(%s)', arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter)));

                            else
                                xlabel(g, sprintf('%s', arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter)));
                            end
                        end
                    end
                    ylabel(g, sprintf('%s [%s]', ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2}));
                end
                
                arSpacedAxisLimits(g);
            end
            countz = 0;
            for jz = iz
                countz = countz + 1;
                g = ar.model(jm).plot(jplot).gz(jz);
                
                set(g,'UserData',...
                    struct('jm',jm,'jplot',jplot,'jx',jx, ...
                    'dLink',ar.model(jm).plot(jplot).dLink, ...
                    'cLink',clinks{jm}{jplot}, ...
                    'model_name',ar.model(jm).name, ...
                    'plot_name',ar.model(jm).plot(jplot).name ...
                    ))
                
                if(~fastPlotTmp)
                    hold(g, 'off');
                    
                    title(g, arNameTrafo(ar.model(jm).z{jz}));
                    
                    if(nu == 0 && nx == 0 && jz == iz(1))
                        if(plot_x_collected)
                            legend(g, cclegendstyles, arNameTrafo([ar.model(jm).u(iu) ar.model(jm).x(ix) ar.model(jm).z(iz)]))
                        else
                            if((~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                                if(~ar.model(jm).plot(jplot).doseresponse)
                                    if(length(ar.model(jm).plot(jplot).dLink)>1)
                                        legend(g, cclegendstyles, arNameTrafo(ar.model(jm).plot(jplot).condition))
                                    end
                                else
                                    legendtmp = {};
                                    ccount = 1;
                                    for jt=1:length(times)
                                        if(~isempty(conditions))
                                            for jc = 1:length(conditions)
                                                legendtmp{ccount} = sprintf('t=%g%s : %s', times(jt), ar.model(jm).tUnits{2}, conditions{jc}); %#ok<AGROW>
                                                ccount = ccount + 1;
                                            end
                                        else
                                            legendtmp{ccount} = sprintf('t=%g%s', times(jt), ar.model(jm).tUnits{2}); %#ok<AGROW>
                                            ccount = ccount + 1;
                                        end
                                    end
                                    legend(g, cclegendstyles, arNameTrafo(legendtmp))
                                end
                            end
                        end
                    end
                    
                    if(countz+nu+nx == (nrows-1)*ncols + 1)
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                        else
                            if(logplotting_xaxis)
                                xlabel(g, sprintf('log_{10}(%s)', arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter)));
                            else
                                xlabel(g, sprintf('%s', arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter))); 
                            end
                        end
                    end
                    ylabel(g, sprintf('%s [%s]', ar.model(jm).zUnits{jz,3}, ar.model(jm).zUnits{jz,2}));
                end
                arSpacedAxisLimits(g);
            end

            if(saveToFile)
                if(ar.config.ploterrors == -1)
                    ar.model(jm).plot(jplot).savePath_FigXCI = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                else
                    ar.model(jm).plot(jplot).savePath_FigX = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                end
            end
            
            figcount = figcount + 1;
        else
            try %#ok<TRYNC>
                close(ar.model(jm).plot(jplot).fighandel_x)
            end
            ar.model(jm).plot(jplot).fighandel_x = [];
        end
    end
end

function [t, u, x, z, ulb, uub, xlb, xub, zlb, zub, jc, dxdt] = getData(jm, jd)
global ar

if(jd~=0)
    jc = ar.model(jm).data(jd).cLink;
else
    jc = 1;
end
t = ar.model(jm).condition(jc).tFine;
u = ar.model(jm).condition(jc).uFineSimu;
x = ar.model(jm).condition(jc).xFineSimu;
z = ar.model(jm).condition(jc).zFineSimu;
if(isfield(ar.model(jm).condition(jc), 'xExpUB'))
    ulb = ar.model(jm).condition(jc).uFineLB;
    uub = ar.model(jm).condition(jc).uFineUB;
    xlb = ar.model(jm).condition(jc).xFineLB;
    xub = ar.model(jm).condition(jc).xFineUB;
    zlb = ar.model(jm).condition(jc).zFineLB;
    zub = ar.model(jm).condition(jc).zFineUB;
else
    ulb = [];
    uub = [];
    xlb = [];
    xub = [];
    zlb = [];
    zub = [];
end
dxdt = ar.model(jm).condition(jc).dxdt;
dxdt(ar.model(jm).condition(jc).qSteadyState==0) = nan;



function [t, u, lb, ub, zero_break] = getDataDoseResponseU(jm, ju, ds, ttime, dLink, logplotting_xaxis)
global ar

zero_break = [];

ccount = 1;
for jd = ds
    for jc = 1:length(ar.model(jm).data(jd).condition)
        if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
            jcondi = jc;
        end
    end
    
    jc = ar.model(jm).data(jd).cLink;
    qt = ar.model(jm).condition(jc).tExp == ttime;
    for jt = find(qt')
        if(logplotting_xaxis)
            t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(jcondi).value)); %#ok<AGROW>
        else
            t(ccount,1) = str2double(ar.model(jm).data(jd).condition(jcondi).value); %#ok<AGROW>
        end
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = dLink
                if(logplotting_xaxis)
                    if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value))))
                        doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value)); %#ok<AGROW>
                    end
                else
                    doses(end+1) = str2double(ar.model(jm).data(jd2).condition(jcondi).value); %#ok<AGROW>
                end
            end
            doses = unique(doses); %R2013a compatible
            if(length(doses)>1)
                t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            else
                t(ccount,1) = doses(1) - 0.1; %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            end
        end
        u(ccount,1) = ar.model(jm).condition(jc).uExpSimu(jt,ju); %#ok<AGROW>
        
        if(isfield(ar.model(jm).condition(jc), 'uExpUB'))
            lb(ccount,1) = ar.model(jm).condition(jc).uExpLB(jt,ju); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).condition(jc).uExpUB(jt,ju); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end
[t,it] = sort(t);
u = u(it);
if(~isempty(lb))
    lb = lb(it);
    ub = ub(it);
end


function [t, x, lb, ub, zero_break] = getDataDoseResponseX(jm, jx, ds, ttime, dLink, logplotting_xaxis)
global ar

zero_break = [];

ccount = 1;
for jd = ds
    for jc = 1:length(ar.model(jm).data(jd).condition)
        if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
            jcondi = jc;
        end
    end
    
    jc = ar.model(jm).data(jd).cLink;
    qt = ar.model(jm).condition(jc).tExp == ttime;
    for jt = find(qt')
        if(logplotting_xaxis)
            t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(jcondi).value)); %#ok<AGROW>
        else
            t(ccount,1) = str2double(ar.model(jm).data(jd).condition(jcondi).value); %#ok<AGROW>
        end
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = dLink
                if(logplotting_xaxis)
                    if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value))))
                        doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value)); %#ok<AGROW>
                    end
                else
                    doses(end+1) = str2double(ar.model(jm).data(jd2).condition(jcondi).value); %#ok<AGROW>
                end
            end
			doses = unique(doses); %R2013a compatible
            if(length(doses)>1)
                t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            else
                t(ccount,1) = doses(1) - 0.1; %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            end
        end
        x(ccount,1) = ar.model(jm).condition(jc).xExpSimu(jt,jx); %#ok<AGROW>
        
        if(isfield(ar.model(jm).condition(jc), 'xExpUB'))
            lb(ccount,1) = ar.model(jm).condition(jc).xExpLB(jt,jx); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).condition(jc).xExpUB(jt,jx); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end

[t,it] = sort(t);
x = x(it);
if(~isempty(lb))
    lb = lb(it);
    ub = ub(it);
end


function [t, z, lb, ub, zero_break] = getDataDoseResponseZ(jm, jz, ds, ttime, dLink, logplotting_xaxis)
global ar

zero_break = [];

ccount = 1;
for jd = ds
    for jc = 1:length(ar.model(jm).data(jd).condition)
        if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
            jcondi = jc;
        end
    end
    
    jc = ar.model(jm).data(jd).cLink;
    qt = ar.model(jm).condition(jc).tExp == ttime;
    for jt = find(qt')
        if(logplotting_xaxis)
            t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(jcondi).value)); %#ok<AGROW>
        else
            t(ccount,1) = str2double(ar.model(jm).data(jd).condition(jcondi).value); %#ok<AGROW>
        end
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = dLink
                if(logplotting_xaxis)
                    if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value))))
                        doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value)); %#ok<AGROW>
                    end
                else
                    doses(end+1) = str2double(ar.model(jm).data(jd2).condition(jcondi).value); %#ok<AGROW>
                end
            end
			doses = unique(doses); %R2013a compatible
            if(length(doses)>1)
                t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            else
                t(ccount,1) = doses(1) - 0.1; %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            end
        end
        z(ccount,1) = ar.model(jm).condition(jc).zExpSimu(jt,jz); %#ok<AGROW>
        
        if(isfield(ar.model(jm).condition(jc), 'zExpUB'))
            lb(ccount,1) = ar.model(jm).condition(jc).zExpLB(jt,jz); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).condition(jc).zExpUB(jt,jz); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end

[t,it] = sort(t);
z = z(it);
if(~isempty(lb))
    lb = lb(it);
    ub = ub(it);
end

%% sub-functions

function C = myLineStyle(n, j)
farben = lines(n);
farben(1,:) = [0 0 0];
C = cell(1,2);
C{1} = 'Color';
C{2} = farben(j,:);



function [h, fastPlotTmp] = myRaiseFigure(m, jplot, figname, figcount, fastPlot)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot(jplot).time = now;
fastPlotTmp = fastPlot;

if(ar.config.ploterrors == -1)
    if(isfield(ar.model(m).plot(jplot), 'fighandel_xCI') && ~isempty(ar.model(m).plot(jplot).fighandel_xCI) && ...
            ar.model(m).plot(jplot).fighandel_xCI ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_xCI==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_xCI, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_xCI;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.4+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_xCI = h;
        fastPlotTmp = false;
    end    
else
    if(isfield(ar.model(m).plot(jplot), 'fighandel_x') && ~isempty(ar.model(m).plot(jplot).fighandel_x) && ...
            ar.model(m).plot(jplot).fighandel_x ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_x==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_x, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_x;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.4+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_x = h;
        fastPlotTmp = false;
    end
end
if(~fastPlot)
    clf
end



function savePath = mySaveFigure(h, name)
global ar
if(ar.config.ploterrors == -1)
    savePath = [arSave '/FiguresCI/X'];
else
    savePath = [arSave '/Figures/X'];
end

if(~exist(savePath, 'dir'))
    mkdir(savePath)
end

if(length(name)>50)
    name = name(1:50);
end

savePath = arPathConvert([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc2', savePath);
if(ispc)
    print('-dpdf', savePath);
elseif(ismac)
    system(['/usr/local/bin/ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
else
    system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
end







function [ncols, nrows, nu, nx, nz, iu, ix, iz] = myColsAndRows(jm)
global ar
if(~isfield(ar.model(jm), 'qPlotU'))
    nu = size(ar.model(jm).u, 2);
    iu = 1:nu;
else
    nu = sum(ar.model(jm).qPlotU);
    iu = find(ar.model(jm).qPlotU);
end
if(~isfield(ar.model(jm), 'qPlotX'))
    nx = size(ar.model(jm).x, 2);
    ix = 1:nx;
else
    nx = sum(ar.model(jm).qPlotX);
    ix = find(ar.model(jm).qPlotX);
end
if(~isfield(ar.model(jm), 'qPlotZ'))
    nz = size(ar.model(jm).z, 2);
    iz = 1:nz;
else
    nz = sum(ar.model(jm).qPlotZ);
    iz = find(ar.model(jm).qPlotZ);
end
[nrows, ncols] = arNtoColsAndRows(nu+nx+nz);




