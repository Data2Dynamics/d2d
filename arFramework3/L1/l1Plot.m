% Plot L1 scan summary

function l1Plot

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    if(~isfield(ar,'L1jks') || isempty(ar.L1jks))
        error('please initialize by l1Init, run l1Scan, l1Unpen, and l1SelectOpt')
    end
end

jks = ar.L1jks;
linv = ar.L1linv;
ps = ar.L1ps;
chi2s_unpen = ar.L1chi2s_unpen;
final_ind = ar.L1final_ind;
parsgt0 = ar.L1parsgt0;
signifmat = ar.L1signifmat;

l = arNameTrafo(ar.pLabel(jks));

figure

subplot(2,2,[2 4])
im = imagesc(ps(:,jks)');
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

[row,col] = find(abs(ps(:,jks)) > 1e-6);
parbar = zeros(1,length(jks));
for i = 1:length(jks)
    myind = col == i;
    if sum(myind) > 0
        parbar(i) = max(row(myind));
    end
end

hb = barh(1:size(im_cdata,1),linvlog(parbar)+.5*ltick,'FaceColor','none','BarWidth',1,'BaseValue',linvlog(1)-.5*ltick);

colmax = min([5 ceil(max(max(abs(ps))))]);
set(gca,'CLim',[-colmax colmax])
set(gca,'YLim',[0 length(jks)+1])
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
set(gca,'YTick',1:length(jks),'YTickLabel',l)

xlabel('log_{10}(\lambda)')
myxlim = get(gca,'XLim');

subplot(2,2,1)
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
plot(linvlogstep+.5*ltick,parsgt0step)
hold on
plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,get(gca,'YLim'),'k:')
xlabel('log_{10}(\lambda)')
ylabel('No. of cell-type specific parameters')
xlim(myxlim)

subplot(2,2,3)
a = semilogy(linvlogstep+.5*ltick,chi2s_unpenstep-chi2s_unpenstep(1)+1);
hold on
set(gca,'YLim',[.4 10^ceil(log10(max(signifmatstep)))])
xlabel('log_{10}(\lambda)')
ylabel('Likelihood ratio')
myylim = get(gca,'YLim');
b = plot(linvlogstep,1-(signifmatstep-chi2s_unpenstep+chi2s_unpenstep(1)),'r--');
c = plot([linvlog(final_ind) linvlog(final_ind)]+.5*ltick,myylim,'k:');
xlim(myxlim)
box on
legend([a b c],{'Test statistic D + 1','Statistical threshold','Parsimonious model'})