% l1Plot([eraseStr])
% 
% Plot L1 scan summary
% 
% eraseStr      An optional string to be erased from parameter names
%               i.e. from ar.pLabel

function l1Plot(eraseStr)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~isfield(ar,'L1jks') || isempty(ar.L1jks))
    error('please initialize by l1Init, run l1Scan, l1Unpen, and l1SelectOpt')
end

jks = ar.L1jks;
linv = ar.linv;
ps = ar.L1ps;
chi2s_unpen = ar.L1chi2s_unpen;
chi2s_lam0 = ar.L1lam0chi2s;
final_ind = ar.L1final_ind;
parsgt0 = ar.L1parsgt0;
signifmat = ar.L1signifmat;

l = arNameTrafo(ar.pLabel(jks));

% Calculate differences of fc parameters
ps2 = []; jks2 = [];
if ar.L1DiffPen_activate
    diffIds = ar.L1DiffPen_diffs;
    eqParMap = ar.L1DiffPen_constrSwitches;
    ps2 = bsxfun(@minus, ps(:,diffIds(:,1)), ps(:,diffIds(:,2)));
    ps2(ar.L1ps_unpen(:,eqParMap) == 1) = 0; % remove diffs that are turned of by iseq from plot
    jks2 = 1:size(diffIds,1);
    ps2Labels = append('Diff_',ar.pLabel(diffIds(:,1)),'_',ar.pLabel(diffIds(:,2)));
    l = [l,arNameTrafo(ps2Labels)]; %%% CUSTOM
end

if (exist('eraseStr','var') && ~isempty(eraseStr) && ischar(eraseStr))
    for i = 1:length(l)
        l{i} = erase(l{i},eraseStr);
    end
end

figure

subplot(2,5,[4 5 9 10])
im = imagesc([ps(:,jks),ps2]');
im_cdata = im.CData;
delete(im)
linvlog = -log10(linv);
ltick = linvlog(2)-linvlog(1);
for j = 1:size(im_cdata,1)
    for i = 1:size(im_cdata,2)
        patch([linvlog(i)-.5*ltick linvlog(i)+.5*ltick linvlog(i)+.5*ltick linvlog(i)-.5*ltick],[j-.5 j-.5 j+.5 j+.5],im_cdata(j,i), 'EdgeColor', 'none')
        hold on
    end
end
xlim([linvlog(1)-.5*ltick linvlog(end)+.5*ltick])
ylim([.5 j+.5])
set(gca,'YDir','reverse')
box on
mypos = get(gca,'Pos');
hold on

% [row,col] = find(abs(ps(:,jks)) > ar.L1thresh);
% parbar = zeros(1,length(jks));
% for i = 1:length(jks)
%     myind = col == i;
%     if sum(myind) > 0
%         parbar(i) = max(row(myind));
%     end
% end

rc = abs([ps(:,jks),ps2]) > ar.L1thresh;
for i = 1:length([jks,jks2])
    begin = NaN;
    for j = 1:(length(linv)-1)
        if rc(j,i) && isnan(begin)
            begin = j;
        elseif ~rc(j,i) && ~isnan(begin)
            stop = j-1;
            patch([linvlog(begin)-.5*ltick linvlog(stop)+.5*ltick linvlog(stop)+.5*ltick linvlog(begin)-.5*ltick],...
                [i-.48 i-.48 i+.48 i+.48],1, 'FaceColor', 'none','EdgeColor','black')
            begin = NaN;
        end
    end
    if ~isnan(begin)
        stop = length(linv);
        patch([linvlog(begin)-.5*ltick linvlog(stop)+.5*ltick linvlog(stop)+.5*ltick linvlog(begin)-.5*ltick],...
                [i-.48 i-.48 i+.48 i+.48],1, 'FaceColor', 'none','EdgeColor','black')
    end
end

x = 1:size(im_cdata,1);

% hb = barh(x(parbar>0),linvlog(parbar(parbar>0))+.5*ltick,...
%     'FaceColor','none','BarWidth',1,'BaseValue',linvlog(1)-.5*ltick);

colmax = min([5 ceil(max(max(abs([ps,ps2]))))]); 
set(gca,'CLim',[-colmax colmax])
set(gca,'YLim',[0 length([jks jks2])+1])
plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,get(gca,'YLim'),'k:')
colmap = [];
for i = 1:250
    colmap = [colmap; hsv2rgb([0 (1-(i-1)/250)^.7 1])];
end
colmap = [colmap; [1 1 1]];
for i = 1:250
    colmap = [colmap; hsv2rgb([2/3 (i/250)^.7 1])];
end
colormap(gca,colmap)
colorbar('Location','EastOutside')
set(gca,'Pos',mypos)

set(gca,'YDir','Reverse')
set(gca,'YTick',1:length([jks jks2]),'YTickLabel',l)

xlabel('log_{10}(\lambda)')
myxlim = get(gca,'XLim');

subplot(2,5,[1 2])
steppars = 1:length(linv)-1;
linvlogstep = linvlog;
chi2s_unpenstep = chi2s_unpen;
signifmatstep = signifmat;
parsgt0step = parsgt0;

for i = 1:length(steppars)
    parsgt0step = parsgt0step([1:steppars(i)-1+i steppars(i)+i steppars(i)+i:end]);
    chi2s_unpenstep = chi2s_unpenstep([1:steppars(i)-1+i steppars(i)+i steppars(i)+i:end]);
    signifmatstep = signifmatstep([1:steppars(i)-1+i steppars(i)+i steppars(i)+i:end]);
    linvlogstep = linvlogstep([1:steppars(i)-1+i steppars(i)+i-1 steppars(i)+i:end]);

end

switch ar.L1subtype(jks(1))
    case 2
        yyaxis right
        plot(linvlog,ar.gamma,'+')
        hold on
        ylabel('Adapt. Lasso Exponent')
        ylim([-0.05, 1.05])
        yyaxis left
    case 3
        yyaxis right
        plot(linvlog,ar.nu,'+')
        hold on
        ylabel('L_q Exponent')
        ylim([-0.05, 1.05])
        yyaxis left
    case 4
        yyaxis right
        plot(linvlog,ar.alpharange,'+')
        hold on
        ylabel('Elasticity')
        ylim([-0.05, 1.05])
        yyaxis left
end

plot(linvlogstep+.5*ltick,parsgt0step)
hold on
plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,get(gca,'YLim'),'k:')
ylabel('No. of cell-type specific parameters')


xlabel('log_{10}(\lambda)')

xlim(myxlim)

subplot(2,5,[6 7])
% maxy = max(max(chi2s_unpenstep-chi2s_lam0+1,1-(signifmatstep-chi2s_unpenstep+chi2s_lam0)));
% miny = min(min(chi2s_unpenstep-chi2s_lam0+1,1-(signifmatstep-chi2s_unpenstep+chi2s_lam0)));
% a = semilogy(linvlogstep+.5*ltick,chi2s_unpenstep-chi2s_unpenstep(1)+1);
% hold on
% 
% set(gca,'YLim',[min(.4,miny) 10^ceil(log10(maxy))])
% 
% xlabel('log_{10}(\lambda)')
% ylabel('Likelihood ratio')
% myylim = get(gca,'YLim');
% b = plot(linvlogstep,1-(signifmatstep-chi2s_unpenstep+chi2s_lam0),'r--');
% c = plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,myylim,'k:');
% xlim(myxlim)
% box on
% legend([a b c],{'Test statistic D + 1','Statistical threshold','Parsimonious model'})

if strcmpi(ar.L1seltype,'LRT')
    maxy = max(max(chi2s_unpenstep-chi2s_lam0+1,1-(signifmatstep-chi2s_unpenstep+chi2s_lam0)));
    a = semilogy(linvlogstep+.5*ltick,chi2s_unpenstep-chi2s_lam0+1);
    hold on
    set(gca,'YLim',[.4 10^ceil(log10(maxy))])
    xlabel('log_{10}(\lambda)')
    ylabel('Likelihood ratio')
    myylim = get(gca,'YLim');
    b = plot(linvlogstep,1-(signifmatstep-chi2s_unpenstep+chi2s_lam0),'r--');
    c = plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,myylim,'k:');
    xlim(myxlim)
    box on
    legend([a b c],{'Test statistic D + 1','Statistical threshold','Parsimonious model'})
    hold off
elseif strcmpi(ar.L1seltype,'BIC')
    maxy = max(signifmatstep);
    miny = min(signifmatstep);
    extension= 0.1 * ((maxy) - (miny));
    b = plot(linvlogstep,signifmatstep,'k-');
    myylim = get(gca,'YLim');
    hold on
    c = plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,myylim,'k:');
    xlim(myxlim)
    set(gca,'YLim',[miny-extension maxy+extension])
    xlabel('log_{10}(\lambda)')
    ylabel('Bayesian Information')
    legend([b c],{'BIC','Parsimonious Model'})
    hold off
end