% pleDebug2(figure_number)
% plots ple parameters

function pleDebug2(figure_number)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end

if(~exist('figure_number', 'var'))
    figure_number = 1;
end

%% Plot Step Statistics

count = 0;
for jk=1:length(ar.ple.ps)
    if(~isempty(ar.ple.ps{jk}))
        figure(figure_number+count)
        
        notjk = 1:length(ar.ple.p);
        notjk = notjk~=jk;
        
        subplot(1,1,1)
        plot(ar.ple.ps{jk}(:,jk), ar.ple.ps{jk}(:,notjk), 'x-')
        hold on
        plot(ar.ple.ps{jk}(:,jk), ar.ple.psinit{jk}(:,notjk), 'x--')
        plot([0 0]+ar.ple.p(jk), ylim, 'k');
        hold off
        ylabel('parameters')
        
        count = count + 1;
    end
end