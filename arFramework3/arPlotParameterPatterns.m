function arPlotParameterPatterns(ps, jks)

global ar

if(exist('ps','var'))
    nbest = size(ps,1);
else
    nbest = 0;
end
if(~exist('jks','var'))
    jks = 1:length(ar.p);
end

hs = [];
hlabel = {};

figure(1); clf;
C = jet(nbest);
C = bsxfun(@rdivide, C, sqrt(sum(C.^2,2)));
patch([ar.lb(jks) fliplr(ar.ub(jks))], [1:length(jks) length(jks):-1:1], ...
    -1*ones(size([1:length(jks) length(jks):-1:1])), ...
    'FaceColor', [.8 .8 .8], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
hold on
for j=1:nbest
    h = plot(ps(j,jks), 1:length(jks), 'Color', C(j,:));
    if(j==nbest)
        hs(end+1) = h; %#ok<AGROW>
        hlabel{end+1} = 'fits'; %#ok<AGROW>
    end
end
h = plot(ar.p(jks), 1:length(jks), 'ko');
hs(end+1) = h;
hlabel{end+1} = 'current value';
if(isfield(ar, 'pTrue'))
    h = plot(ar.pTrue(jks), 1:length(jks), 'k--');
    hs(end+1) = h;
    hlabel{end+1} = 'true value';
end
hold off
title(sprintf('parameter differences between %i best fits', nbest));
xlabel('parameter values')
set(gca, 'YTick', 1:length(jks));
set(gca, 'YTickLabel', cellfun(@(x,y) sprintf('%s #%i',x,y), ...
    arNameTrafo(ar.pLabel(jks)), num2cell(jks), 'UniformOutput',false));
set(gca, 'YDir','reverse');
ylim([0 length(jks)+1])
xrange = 0.01*(max(ar.ub(jks)) - min(ar.lb(jks)));
xlim([min(ar.lb(jks))-xrange max(ar.ub(jks))+xrange])
grid on
legend(hs, hlabel,'Location','Best');
