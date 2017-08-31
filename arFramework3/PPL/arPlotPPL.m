% plot prediction profile likelihood
%
% arPlotPPL

function arPlotPPL(m, c, ix, t, takeY, subs_para)

global ar
qLog10 = ar.ppl.qLog10;

if(~exist('subs_para','var'))
    subs_para = true;
end

if(~exist('takeY','var') || isempty(takeY))
   takeY = true; 
end
if(takeY)
    data_cond = 'data';
else
    data_cond = 'condition';
end

if(length(ix)>1)
    error('length of ix must be = 1');
end
if(length(t)>1)
    error('length of t must be = 1');
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

% qx = ar.model(m).(data_cond)(c).ppl.ix == ix;
qt = ar.model(m).(data_cond)(c).ppl.tstart(:,ix) == t;

xtrial = squeeze(ar.model(m).(data_cond)(c).ppl.xtrial(qt,ix,:));
xfit = squeeze(ar.model(m).(data_cond)(c).ppl.xfit(qt,ix,:));
ppl = squeeze(ar.model(m).(data_cond)(c).ppl.vpl(qt,ix,:));
if(ar.ppl.options.doPPL || sum(isnan(ppl)) == length(ppl))
    ppl = squeeze(ar.model(m).(data_cond)(c).ppl.ppl(qt,ix,:));
end
ps = squeeze(ar.model(m).(data_cond)(c).ppl.ps(qt,ix,:,:));

if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1)
    ppl = 2*ar.ndata*log(sqrt(2*pi)) + ppl;
elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
    ppl = 2*ar.ndata_err*log(sqrt(2*pi)) + ppl;
end

farben = lines(length(ar.pLabel));
zeichen = {'-', '--', '-.'};
ps_label_count = 5;
ps_var_min_label = 1e-2;

figure(1);

subplot(2,1,1)

plot(xfit, ppl, 'k')
hold on
plot(xtrial, ppl, 'k--')
plot(xtrial(ceil(length(xtrial)/2)), min(ppl), '*', 'Color', [.5 .5 .5], 'LineWidth', 1, 'MarkerSize', 8)
plot(xlim, [0 0]+min(ppl)+ar.ppl.dchi2, 'r--')
text(mean(xlim), min(ppl)+ar.ppl.dchi2, sprintf('%2i%% (ndof = %i)', (1-ar.ppl.alpha_level)*100, ...
    ar.ppl.ndof), 'Color', 'r', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
hold off
legend('PPL','VPL');
title(sprintf('%s (m=%i, c=%i, t=%g)', arNameTrafo(ar.model(m).x{ix}), m, c, t));
xlimtmp = xlim;

if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1) || ...
        (ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0) )
    ylabel('-2*log(PL)');
else
    ylabel('\chi^2_{PL}');
end

g = subplot(2,1,2);

qFit = ar.qFit == 1;
legendstmp = cell(1,sum(qFit));
legendstmplines = zeros(1,sum(qFit));
ps_var = zeros(1,sum(qFit));

ccount = 1;
for j=find(qFit)
    qisnonan = ~isnan(xfit);
    ps_var(ccount) = max(ps(qisnonan,j)) - min(ps(qisnonan,j));
    zeichenindex = mod(floor((j-1)/7)+1, 3)+1;
    
    if(subs_para)
        medianp = ar.p(j);
    else
        medianp = 0;
    end
    
    line_s = plot(xfit(qisnonan), ps(qisnonan,j)-medianp, [zeichen{zeichenindex}], 'color', farben(j,:));
    hold on
    
    legendstmp{ccount} = arNameTrafo(ar.pLabel{j}); %#ok<*agrow>
    legendstmplines(ccount) = line_s;
    ccount = ccount + 1;
end
hold off
if(~isempty(legendstmplines))
    if(ps_label_count<length(legendstmp))
        [ps_varsorted, i_largest_std] = sort(ps_var, 2, 'descend');
        i_largest_std = i_largest_std(1:ps_label_count);
        qi_largest_std = ps_varsorted(1:ps_label_count) > ps_var_min_label;
        i_largest_std = i_largest_std(qi_largest_std);
        if(~isempty(i_largest_std))
            legend(legendstmplines(i_largest_std), legendstmp(i_largest_std), 'location', 'best')
        end
    else
        legend(legendstmplines, legendstmp, 'location', 'best')
    end
end
xlim(xlimtmp);
ylabel({'change of parameters'})
if(qLog10)
    xlabel(['log_{10}(' arNameTrafo(ar.model(m).x{ix}) ')'])
else
    xlabel(arNameTrafo(ar.model(m).x{ix}))
end
grid(g,'on');

