% pleDebug(figure_number)

function pleDebug(figure_number)
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

        subplot(4,2,1)
%         tmpchi2 = diff(pleGlobals.chi2s{jk});
%         tmpchi2 = [tmpchi2(1:floor(length(tmpchi2)/2)) 0 tmpchi2((floor(length(tmpchi2)/2)+1):end)];
        tmpchi2 =  pleGlobals.chi2sinit{jk}(1:floor(length(pleGlobals.chi2s{jk})/2)) - pleGlobals.chi2s{jk}(2:ceil(length(pleGlobals.chi2s{jk})/2));
        tmpchi2 = [tmpchi2 0 pleGlobals.chi2sinit{jk}((ceil(length(pleGlobals.chi2s{jk})/2)+1):end) - pleGlobals.chi2s{jk}(ceil(length(pleGlobals.chi2s{jk})/2):(end-1))];
        plot(pleGlobals.ps{jk}(:,jk), tmpchi2, '*-')
        hold on
        plot(xlim, [0 0]+pleGlobals.relchi2stepincrease(jk)*pleGlobals.dchi2, 'r--')
%         plot(xlim, [0 0]-pleGlobals.relchi2stepincrease(jk)*pleGlobals.dchi2, 'r--')
        plot([0 0]+pleGlobals.p(jk), ylim, 'k');
        hold off
        ylabel('\Delta\chi^2(step)')
        
        subplot(4,2,3)
        plot(pleGlobals.ps{jk}(:,jk), pleGlobals.chi2sinit{jk} - ...
            pleGlobals.chi2s{jk}, '*-')
        hold on
        q_fiterror = pleGlobals.chi2sinit{jk}-pleGlobals.chi2s{jk} < 0;
        plot(pleGlobals.ps{jk}(q_fiterror,jk), pleGlobals.chi2sinit{jk}(q_fiterror) - ...
            pleGlobals.chi2s{jk}(q_fiterror), 'r*')        
        plot([0 0]+pleGlobals.p(jk), ylim, 'k');
        hold off
        ylabel('\chi^2(init-fit)')
        
        subplot(4,2,5)
        semilogy(pleGlobals.ps{jk}(:,jk), sqrt(sum(pleGlobals.psinitstep{jk}.^2,2)), '*-')
        hold on
        semilogy(pleGlobals.ps{jk}(:,jk), abs(pleGlobals.psinitstep{jk}(:,jk)), 'rx-')
        semilogy(xlim, [0 0] + pleGlobals.maxstepsize(jk), 'r--')
        semilogy(xlim, [0 0] + pleGlobals.minstepsize(jk), 'r--')
        plot([0 0]+pleGlobals.p(jk), ylim, 'k');
        hold off
        ylabel('step size')
        
%         subplot(4,2,7)

%         subplot(4,2,[2 4])        
%         plot(pleGlobals.ps{jk}(:,jk), pleGlobals.betas{jk})
%         hold on
%         plot([0 0]+pleGlobals.p(jk), ylim, 'k');
%         hold off
%         ylabel('beta')
% 
%         subplot(4,2,[6 8])
%         plot(pleGlobals.ps{jk}(:,jk), pleGlobals.alphas{jk})
%         hold on
%         plot([0 0]+pleGlobals.p(jk), ylim, 'k');
%         hold off
%         ylabel('alpha')
        
        count = count + 1;
    end
end