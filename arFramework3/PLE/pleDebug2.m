% pleDebug2(figure_number)

function pleDebug2(figure_number)
global pleGlobals;

if(isempty(pleGlobals))
    error('perform ple before usage');
end
if(isempty(pleGlobals.ps))
    return
end

if(~exist('figure_number', 'var'))
    figure_number = 1;
end

%% Plot Step Statistics

count = 0;
for jk=1:length(pleGlobals.ps)
    if(~isempty(pleGlobals.ps{jk}))
        figure(figure_number+count)
        
        notjk = 1:length(pleGlobals.p);
        notjk = notjk~=jk;
        
        subplot(1,1,1)
        plot(pleGlobals.ps{jk}(:,jk), pleGlobals.ps{jk}(:,notjk), 'x-')
        hold on
        plot(pleGlobals.ps{jk}(:,jk), pleGlobals.psinit{jk}(:,notjk), 'x--')
        plot([0 0]+pleGlobals.p(jk), ylim, 'k');
        hold off
        ylabel('parameters')
        
        count = count + 1;
    end
end