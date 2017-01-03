% arPlotClockTimes
% 
%   Show execution times of threads and conditions.

function arPlotClockTimes

global ar

figure(1)

c = 1e-3;

plot([0 ar.stop]*c, [0 0], '*-k');
hold on

labels = {'overall'};

for jt=1:length(ar.config.threads)
    for jmc=1:length(ar.config.threads(jt).ms)
        C = arLineMarkersAndColors(jmc, length(ar.config.threads(jt).ms), [], '*','-');
        jm = ar.config.threads(jt).ms(jmc)+1;
        jc = ar.config.threads(jt).cs(jmc)+1;
        plot([ar.model(jm).condition(jc).start ...
            ar.model(jm).condition(jc).stop_data ...
            ar.model(jm).condition(jc).stop]*c, ...
            [jt jt jt], C{:});
        hold on
    end
    labels{jt+1} = sprintf('thread %i (#mc=%i, #d=%i)', jt, ...
            ar.config.threads(jt).n, ar.config.threads(jt).nd);
end
hold off
ylim([-0.2 ar.config.nParallel+0.2])
set(gca, 'yTick', 0:ar.config.nParallel);
set(gca, 'yTickLabel', labels);
xlabel('time [milli seconds]');
grid on