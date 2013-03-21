function arPlotClockTimes

global ar

figure(1)

ccount = 0;

c = 1e-3;

plot([0 ar.stop]*c, [ccount ccount], '*-k');
hold on

labels = {'overall'};

for jm=1:length(ar.model)
    for jc=1:length(ar.model(jm).condition)
        ccount = ccount + 1;
        plot([ar.model(jm).condition(jc).start ...
            ar.model(jm).condition(jc).stop_data ...
            ar.model(jm).condition(jc).stop]*c, ...
            [ccount ccount ccount], '*-b');
        labels{ccount+1} = sprintf('m=%i c=%i (#d=%i)', jm, jc, ...
            length(ar.model(jm).condition(jc).dLink));
    end
end
hold off
ylim([-0.5 ccount+0.5])
set(gca, 'yTick', 0:ccount);
set(gca, 'yTickLabel', labels);
xlabel('time [milli seconds]');
grid on