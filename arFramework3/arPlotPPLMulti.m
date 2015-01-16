% plot prediction profile likelihood
%
% arPlotPPL

function arPlotPPLMulti(m, c, vpl, filled, dosurf)

global ar

if(~exist('vpl','var'))
    vpl = false;
end
if(~exist('filled','var'))
    filled = false;
end
if(~exist('dosurf','var'))
    dosurf = false;
end

qLog10 = ar.ppl.qLog10;
[nrows, ncols] = arNtoColsAndRows(length(ar.model(m).condition(c).ppl.ix));

chi2 = ar.chi2fit;
if(ar.config.fiterrors == 1)
    chi2 = 2*ar.ndata*log(sqrt(2*pi)) + chi2;
end

arSimu(false, true);

figure(1);
clf;

for jx=1:length(ar.model(m).condition(c).ppl.ix)
    ix = ar.model(m).condition(c).ppl.ix(jx);
    g = subplot(nrows, ncols, jx);
    arSubplotStyle(g)
    
    t = ar.model(m).condition(c).ppl.t;
    if(vpl)
        x = squeeze(ar.model(m).condition(c).ppl.xtrial(:,jx,:));
    else
        x = squeeze(ar.model(m).condition(c).ppl.xfit(:,jx,:));
    end
    ppl = squeeze(ar.model(m).condition(c).ppl.ppl(:,jx,:));
    
    if(ar.config.fiterrors == 1)
        ppl = 2*ar.ndata*log(sqrt(2*pi)) + ppl;
    end
    
    if(qLog10)
        x = 10.^x;
    end
    
    tFine = ar.model(m).condition(c).tFine(ar.model(m).condition(c).tFine>=min(t) & ...
        ar.model(m).condition(c).tFine<=max(t));
    xFine = linspace(min(x(:)), max(x(:)), ar.config.nFinePoints);
    [tgrid, xgrid] = meshgrid(tFine, xFine);
    pplgrid = nan(size(tgrid));
    
    if(~dosurf || filled)
        method = 'v5cubic';
        for jt=1:length(tFine)
            tdiff = t - tFine(jt);
            qback = tdiff<=0;
            qfor = tdiff>0;
            
            iback = find(qback,1,'last');
            ifor = find(qfor,1,'first');
            [~, ifine] = min(abs(ar.model(m).condition(c).tFine - tFine(jt)));
            
            qnonnanback = ~isnan(ppl(iback,:));
            qnonnanfor = ~isnan(ppl(ifor,:));
            
            if(sum(qnonnanfor)>0)
                weigth = 1-abs([tdiff(iback) tdiff(ifor)])/sum(abs([tdiff(iback) tdiff(ifor)]));
                
                xback = x(iback,ceil(length(x)/2));
                xfor = x(ifor,ceil(length(x)/2));
                shift = [ar.model(m).condition(c).xFineSimu(ifine,ix)-xback ...
                    ar.model(m).condition(c).xFineSimu(ifine,ix)-xfor];
                %             shift = [0 0];
                
                if(sum(qnonnanback)>1)
                    pplback = interp1(x(iback,qnonnanback)+shift(1), ppl(iback,qnonnanback), xFine, method);
                elseif(sum(qnonnanback)==1)
                    pplback = 0;
                    weigth(2) = 1;
                end
                if(sum(qnonnanfor)>1)
                    pplfor = interp1(x(ifor,qnonnanfor)+shift(2), ppl(ifor,qnonnanfor), xFine, method);
                elseif(sum(qnonnanfor)==1)
                    pplfor = 0;
                    weigth(1) = 1;
                end
                pplgrid(:,jt) = weigth(1)*pplback + weigth(2)*pplfor;
                
                %             plot(xgrid(:,jt), pplgrid(:,jt))
                %             hold on
                %             plot(xgrid(:,jt), pplback, '--')
                %             plot(xgrid(:,jt), pplfor, '--')
                %             plot(ar.model(m).condition(c).xFineSimu(ifine,ix), min(pplgrid(:,jt)), '*')
                %             title(sprintf('t=%g', tFine(jt)));
                %             hold off
                %             waitforbuttonpress
            else
                pplgrid(:,jt) = interp1(x(iback,qnonnanback), ppl(iback,qnonnanback), xFine, method);
            end
        end
    end
    
    if(~dosurf)
        Cvals = linspace(chi2, chi2+ar.ppl.dchi2*4, 20);
        Cvals = Cvals - (Cvals(2)-Cvals(1))/2;
        
        if(filled)
            map = 1-((1-gray)/2);
            colormap(map);
            contourf(tgrid, xgrid, pplgrid, Cvals,'EdgeColor','none');
            h = colorbar;
            hold(h,'on')
            plot(h, [-1 2], [0 0]+chi2+ar.ppl.dchi2, 'r')
            hold(h,'off')
            hold on
            contour(tgrid, xgrid, pplgrid, chi2+ar.ppl.dchi2, 'r')
            plot(ar.model(m).condition(c).tFine, ar.model(m).condition(c).xFineSimu(:,ix), 'k');
            hold off
            if(~vpl)
                legend('PPL', 'PCI', 'best fit');
            else
                legend('VPL', 'VCI', 'best fit');
            end
        else
            map = gray;
            colormap(map);
            plot(ar.model(m).condition(c).tFine, ar.model(m).condition(c).xFineSimu(:,ix), 'k--');
            hold on
            contour(tgrid, xgrid, pplgrid, Cvals);
            h = colorbar;
            hold(h,'on')
            plot(h, [-1 2], [0 0]+chi2+ar.ppl.dchi2, 'r')
            hold(h,'off')
            contour(tgrid, xgrid, pplgrid, chi2+ar.ppl.dchi2, 'r')
            hold off
            if(~vpl)
                legend('best fit', 'PPL', 'PCI');
            else
                legend('best fit', 'VPL', 'VCI');
            end
        end
    else
        if(filled)
            xfitfine = ar.model(m).condition(c).xFineSimu(:,ix);
            plot3(ar.model(m).condition(c).tFine, xfitfine, ...
                chi2*ones(size(ar.model(m).condition(c).tFine))+0.1, 'k');
            hold on
            
            n = 100;
            [ttmp, xtmp] = meshgrid(linspace(min(tFine), max(tFine), n), ...
                linspace(min(xFine), max(xFine), n));
            ppltmp = interp2(tgrid, xgrid, pplgrid, ttmp, xtmp);
            surf(ttmp, xtmp, ppltmp);
            colormap(cool);
            caxis([chi2 chi2+ar.ppl.dchi2*2]);
            zlim([chi2 chi2+ar.ppl.dchi2*2]);
            h = colorbar;
            hold(h,'on')
            plot(h, [-1 2], [0 0]+chi2+ar.ppl.dchi2, 'r')
            hold(h,'off')
            contour3(tgrid, xgrid, pplgrid, [chi2+ar.ppl.dchi2 chi2+ar.ppl.dchi2], 'r')
            hold off
            if(~vpl)
                legend('best fit', 'PPL', 'PCI');
            else
                legend('best fit', 'VPL', 'VCI');
            end
            
            view(-10, 80);
            camlight('headlight');
            box on
            grid on
        else
            xfitfine = ar.model(m).condition(c).xFineSimu(:,ix);
            h(1) = plot3(ar.model(m).condition(c).tFine, xfitfine, ...
                zeros(size(ar.model(m).condition(c).tFine))+1e-8, 'k--');
            hold on
            
%             n = 100;
%             [ttmp, xtmp] = meshgrid(linspace(min(tFine), max(tFine), n), ...
%                 linspace(min(xFine), max(xFine), n));
%             ppltmp = interp2(tgrid, xgrid, pplgrid, ttmp, xtmp);
%             [~, htmp] = contour3(tgrid, xgrid, pplgrid, [chi2+ar.ppl.dchi2 chi2+ar.ppl.dchi2], 'r');
%             if(~isempty(htmp))
%                 h(3) = htmp(1);
%             end
%             htmp = plot3(ttmp, xtmp, ppltmp,'k');
            for jt = 1:length(t)
                htmp = plot3(t(jt)*ones(size(x(jt,:))), x(jt,:), ppl(jt,:)-chi2,'k');
            end
            h(2) = htmp(1);
            zlim([0 ar.ppl.dchi2]);
            hold off
            if(~vpl)
                legend(h, {'best fit', 'PCI', 'PPL'});
            else
                legend(h, {'best fit', 'VCI', 'VPL'});
            end
            
            view(-10, 80);
            camlight('headlight');
%             box on
            grid on
        end
    end

    title(sprintf('%s', arNameTrafo(ar.model(m).x{ix})));
    
    xlabel(g, sprintf('%s [%s]', ar.model(m).tUnits{3}, ar.model(m).tUnits{2}));
    ylabel(sprintf('%s [%s]', ar.model(m).xUnits{jx,3}, ar.model(m).xUnits{jx,2}))
end


