% res = arPlotSloppiness(whichplot)
% 
%   
%   whichplot   1  X, X', H=X'X
%               2 (default) EV Spectrum and CI from PL
% 
%   
%   xtype       'sy' default ar.model.data.syExpSimu
%               'sres': ar.sres
%   
%   res         result struct

function res = arPlotSloppiness(whichplot,xtype)
if(~exist('whichplot','var') || isempty(whichplot))
    whichplot = 2;
end
if(~exist('xtype','var') || isempty(xtype))
    xtype = 'sy';
end

global ar
global pleGlobals

X = arGetDesignMatrix(xtype);
dop = find(ar.qFit==1);
unlog = find(ar.qLog10~=1);
unlog2 = intersect(dop,unlog);
if ~isempty(unlog2)
    warning('The following parameters are fitted by qLog10~=1: %s',sprintf('%i ',unlog2));
end


% X = X./(sum(abs(X),2)*ones(1,size(X,2)));

H = X(:,dop)'*X(:,dop);
% H = H./(sum(abs(H),2)*ones(1,size(H,2)));

[dummy,rf] = sort(sum(log10(abs(X)),2));


if(sum(whichplot==1)>0)
    figure
    subplot(2,2,2)
    imagesc(log10(abs(X(rf,dop)')))
    title('log_{10}(abs(X(rf,:)''))')
    subplot(2,2,3)
    imagesc(log10(abs(X(rf,dop))))
    title('log_{10}(abs(X(rf,:)))')
    subplot(2,2,4)
    imagesc(log10(abs(H)))
    title('log_{10}(abs(H))')
end

[V,D] = eig(H);
eigen = diag(D);
if(max(abs(imag(eigen)))>1e-8)
    warning('Complex eigenvalues!')
end
eigen = real(eigen);


if(~isempty(pleGlobals))
    ilog = find(pleGlobals.q_log10(dop)==1);
    iunlog = find(pleGlobals.q_log10(dop)==0);
    
    ci_size(ilog) = pleGlobals.conf_ub_point(dop(ilog))-pleGlobals.conf_lb_point(dop(ilog));
    ci_size(iunlog) = log10(pleGlobals.conf_ub_point(dop(iunlog))) - log10(pleGlobals.conf_lb_point(dop(iunlog)));
else
    ci_size = NaN(size(ar.pLabel));
end
nonId = isinf(ci_size);


if(sum(whichplot==2)>0)
    med_eigen = median(log10(eigen(~isinf(eigen))));
    med_ci = median(log10(ci_size(~isinf(ci_size) & ~isnan(ci_size))));
    
    figure;
    set(gca,'YScale','log');
    hold on
    hp1=patch([0.8,2.2,2.2,0.8],10.^(med_eigen+[-3,-3,3,3]),.8*ones(1,3),'EdgeColor','none');
    hp2=patch([2.8,4.2,4.2,2.8],10.^(med_ci+[-.5,-.5,.5,.5]),[0,1,1],'EdgeColor','none');
    set(hp2,'FaceAlpha',.1);

    h1=semilogy(ones(2,1)*eigen','k','LineWidth',2);
    h2=semilogy([3*ones(1,length(ci_size));4*ones(1,length(ci_size))],ones(2,1)*ci_size,'b','LineWidth',2);
    
    titstr = '';
    
    deig = nanmax(log10(abs(eigen)))-nanmin(log10(abs(eigen)));
    titstr = [titstr,sprintf(', EW of H spread over %.2f orders',deig)];    

    CIda = 1; % max. xlim
    if(isempty(pleGlobals))
        CIda = 0;
        titstr = [titstr,', CIs are not available'];
        warning('PLEs are not available, run ple first');
    elseif sum(isnan(ci_size)>0)
        CIda = 0;
        titstr = [titstr,', some CIs are not available'];
        warning('PLEs are not comprehensively available, compute all profiles first.');
    else
        dci = nanmax(log10(ci_size(~isinf(ci_size)))) -nanmin(log10(ci_size(~isinf(ci_size))));
        titstr = [titstr,sprintf(', CI spread over %.2f orders',dci)];
    end
    
    if(sum(nonId)>0)
        titstr = [titstr,sprintf(', %i non id. param (!)',sum(nonId))];
    end
    
    axis tight
    if CIda        
        xlim([0.5,4.5])
    else
        xlim([0.5,2.5])
    end
    yl = log10(get(gca,'Ylim'));
    dyl = yl(2)-yl(1);
    yl2 = yl + [-dyl*.1,dyl*.1];
    set(gca,'YLim',10.^yl2,'XTick',[],'YGrid','on');

    if length(titstr)>60
        fs = 8;
    else
        fs = 10;
    end
    if CIda        
        set(legend([h1(1),hp1,h2(1),hp2],'eigenvalues','6 orders of magn.','CI size','1 order of magn.','Location','NorthEastOutside'),'FontSize',fs)
    else
        set(legend([h1(1),hp1],'eigenvalues','6 orders of magn.','Location','NorthEastOutside'),'FontSize',fs)
    end
    
    if(~isempty(titstr))
        if strcmp(titstr(1),',')==1
            titstr = titstr(2:end);
        end
        if length(titstr)>60
            fs = 8;
        else
            fs = 10;
        end
        title(titstr,'FontSize',fs)
    end
end

if(nargout>0)
    res.dop = dop;
    res.nonId = nonId;
    res.eigen = eigen;
    res.pLabel = ar.pLabel(dop);
    res.H = H;
    res.X = X;
    res.ev = V;
    
    varargout{1} = res;
end
