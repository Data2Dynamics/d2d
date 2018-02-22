% plot prediction profile likelihood
%
% arPlotPPLMulti
% 
% Set likelihood threshold with ar.ppl.options.alpha_level

function arPlotPPLMulti(m, c, takeY, vpl, filled, dosurf)

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
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end
if(~exist('takeY','var') || isempty(takeY))
   takeY = false; 
end
if(takeY)
    data_cond = 'data';
    x_y = 'y';
else
    data_cond = 'condition';
    x_y = 'x';
end
qLog10 = ar.ppl.qLog10;

% optimizer settings (set only once)
fittederrors=ar.config.fiterrors;
ar.config.fiterrors=0;
fit_bkp = ar.qFit(ar.qError==1);
ar.qFit(ar.qError==1)=2;

if(~isfield(ar.model(m).(data_cond)(c),'ppl'))
    error('There is no ppl calculated!\n');
end

xs_withPPL = find(all(~isnan(ar.model(m).(data_cond)(c).ppl.ts_profile),1));

[nrows, ncols] = arNtoColsAndRows(length(xs_withPPL));
arCalcMerit
arSimu(false,true,true)
chi2 = arGetMerit('chi2fit');
if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1)
    chi2 = 2*ar.ndata*log(sqrt(2*pi)) + chi2;
elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
    chi2 = 2*ar.ndata_err*log(sqrt(2*pi)) + chi2;
end

dchi2 = chi2inv(1-ar.ppl.options.alpha_level, 1);

figure(1);
clf;
jx = 0;
for ix=xs_withPPL
    jx = jx+1;
    g = subplot(nrows, ncols, jx);
    arSubplotStyle(g)
    
    t = ar.model(m).(data_cond)(c).ppl.ts_profile(:,ix);
    if(sum(isnan(t)) == length(t))
        continue;
    end
    if(vpl)
        x = squeeze(ar.model(m).(data_cond)(c).ppl.xtrial_profile(:,ix,:));
    else
        x = squeeze(ar.model(m).(data_cond)(c).ppl.xfit_profile(:,ix,:));
    end
    ppl = squeeze(ar.model(m).(data_cond)(c).ppl.ppl_likelihood_profile(:,ix,:));
    if(vpl)
        ppl = squeeze(ar.model(m).(data_cond)(c).ppl.vpl_likelihood_profile(:,ix,:));
    end
    if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1)
        ppl = 2*ar.ndata*log(sqrt(2*pi)) + ppl;
    elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
        ppl = 2*ar.ndata_err*log(sqrt(2*pi)) + ppl;
    end
    
    if(qLog10)
        x = 10.^x;
    end
    
    tFine = ar.model(m).(data_cond)(c).tFine(ar.model(m).(data_cond)(c).tFine>=min(t) & ...
        ar.model(m).(data_cond)(c).tFine<=max(t));
    %tFine = ar.model(m).(data_cond)(c).tFine;
    xFine = linspace(min(x(:)), max(x(:)), ar.config.nFinePoints);
    [tgrid, xgrid] = meshgrid(tFine, xFine);
    pplgrid = nan(size(tgrid));
    
     %if(~dosurf || filled)
        method = 'spline';
        for jt=1:length(tFine)
            tdiff = t - tFine(jt);
            qback = tdiff<=0;
            qfor = tdiff>0;
            
            iback = find(qback,1,'last');
            ifor = find(qfor,1,'first');
            [~, ifine] = min(abs(ar.model(m).(data_cond)(c).tFine - tFine(jt)));
            
            qnonnanback = ~isnan(ppl(iback,:));
            qnonnanfor = ~isnan(ppl(ifor,:));
            
            if(sum(qnonnanfor)>0)
                weigth = 1-abs([tdiff(iback) tdiff(ifor)])/sum(abs([tdiff(iback) tdiff(ifor)]));
                
                xback = x(iback,ceil(length(x)/2));
                xfor = x(ifor,ceil(length(x)/2));
                shift = [ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(ifine,ix)-xback ...
                    ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(ifine,ix)-xfor];
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
                %             plot(ar.model(m).(data_cond)(c).xFineSimu(ifine,ix), min(pplgrid(:,jt)), '*')
                %             title(sprintf('t=%g', tFine(jt)));
                %             hold off
                %             waitforbuttonpress
            else
                pplgrid(:,jt) = interp1(x(iback,qnonnanback), ppl(iback,qnonnanback), xFine, method);
            end
        end
     %end
    
    if(~dosurf)
        Cvals = linspace(chi2, chi2+dchi2*4, 20);
        Cvals = Cvals - (Cvals(2)-Cvals(1))/2;
        
        if(filled)
            map = 1-((1-gray)/2);
            colormap(map);
            contourf(tgrid, xgrid, pplgrid, Cvals,'EdgeColor','none');
            h = colorbar;
%             hold(h,'on')
%             plot(h, [-1 2], [0 0]+chi2+dchi2, 'r')
%             hold(h,'off')
            hold on
            contour(tgrid, xgrid, pplgrid, 'r')%, chi2+dchi2
            plot(ar.model(m).(data_cond)(c).tFine, ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(:,ix), 'k');
            hold off
            if(~vpl)
                legend('PPL', 'PCI', 'best fit');
            else
                legend('VPL', 'VCI', 'best fit');
            end
        else
            map = gray;
            colormap(map);
            plot(ar.model(m).(data_cond)(c).tFine, ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(:,ix), 'k--');
            hold on
            contour(tgrid, xgrid, pplgrid, Cvals);
%             h = colorbar;
%             hold(h,'on')
%             plot(h, [-1 2], [0 0]+chi2+dchi2, 'r')
%             hold(h,'off')
            contour(tgrid, xgrid, pplgrid,'r')%chi2+dchi2,
            hold off
            if(~vpl)
                legend('best fit', 'PPL', 'PCI');
            else
                legend('best fit', 'VPL', 'VCI');
            end
        end
    else
        if(filled)
            xfitfine = ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(:,ix);
            plot3(ar.model(m).(data_cond)(c).tFine, xfitfine, ...
                chi2*ones(size(ar.model(m).(data_cond)(c).tFine))+0.1, 'k');
            hold on
            
            n = 100;
            [ttmp, xtmp] = meshgrid(linspace(min(tFine), max(tFine), n), ...
                linspace(min(xFine), max(xFine), n));
            ppltmp = interp2(tgrid, xgrid, pplgrid, ttmp, xtmp);
            surf(ttmp, xtmp, ppltmp);
            colormap(cool);
            caxis([chi2 chi2+dchi2*2]);
            zlim([chi2 chi2+dchi2*2]);
            h = colorbar;
            contour3(tgrid, xgrid, pplgrid, [chi2+dchi2 chi2+dchi2], 'r')
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
            xfitfine = ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(:,ix);
            h(1) = plot3(ar.model(m).(data_cond)(c).tFine, xfitfine, ...
                zeros(size(ar.model(m).(data_cond)(c).tFine))+1e-8, 'k--');
            hold on
            
            n = 100;
            [ttmp, xtmp] = meshgrid(linspace(min(tFine), max(tFine), n), ...
                linspace(min(xFine), max(xFine), n));
            ppltmp = interp2(tgrid, xgrid, pplgrid-chi2, ttmp, xtmp);
            %[~, htmp] = contour3(tgrid, xgrid, pplgrid-chi2, 'r');%, [chi2+dchi2 chi2+dchi2]
            [~, htmp] = contour3(tgrid, xgrid, pplgrid-chi2, [dchi2 dchi2], 'r');
            if(~isempty(htmp))
                h(3) = htmp(1);
            end
            %htmp = plot3(ttmp, xtmp, ppltmp,'k');
            for jt = 1:length(t)
                htmp = plot3(t(jt)*ones(size(x(jt,:))), x(jt,:), ppl(jt,:)-chi2,'k');
            end
            h(2) = htmp(1);
            zlim([0 dchi2+1]);
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
    set(gcf,'Color','w')
    xlabel(g, sprintf('%s [%s]', ar.model(m).tUnits{3}, ar.model(m).tUnits{2}));
    ylabel(sprintf('%s [%s]', ar.model(m).xUnits{ix,3}, ar.model(m).xUnits{ix,2}))
end

ar.config.fiterrors=fittederrors;
ar.qFit(ar.qError==1)=fit_bkp;
