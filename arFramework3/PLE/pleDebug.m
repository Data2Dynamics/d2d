% pleDebug(figure_number)

function pleDebug(figure_number)

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

        subplot(4,2,1)
%         tmpchi2 = diff(ar.ple.chi2s{jk});
%         tmpchi2 = [tmpchi2(1:floor(length(tmpchi2)/2)) 0 tmpchi2((floor(length(tmpchi2)/2)+1):end)];
        tmpchi2 =  ar.ple.chi2sinit{jk}(1:floor(length(ar.ple.chi2s{jk})/2)) - ar.ple.chi2s{jk}(2:ceil(length(ar.ple.chi2s{jk})/2));
        tmpchi2 = [tmpchi2 0 ar.ple.chi2sinit{jk}((ceil(length(ar.ple.chi2s{jk})/2)+1):end) - ar.ple.chi2s{jk}(ceil(length(ar.ple.chi2s{jk})/2):(end-1))];
        plot(ar.ple.ps{jk}(:,jk), tmpchi2, '*-')
        hold on
        plot(xlim, [0 0]+ar.ple.relchi2stepincrease(jk)*ar.ple.dchi2, 'r--')
%         plot(xlim, [0 0]-ar.ple.relchi2stepincrease(jk)*ar.ple.dchi2, 'r--')
        plot([0 0]+ar.ple.p(jk), ylim, 'k');
        hold off
        ylabel('\Delta\chi^2(step)')
        
        subplot(4,2,3)
        plot(ar.ple.ps{jk}(:,jk), ar.ple.chi2sinit{jk} - ...
            ar.ple.chi2s{jk}, '*-')
        hold on
        q_fiterror = ar.ple.chi2sinit{jk}-ar.ple.chi2s{jk} < 0;
        plot(ar.ple.ps{jk}(q_fiterror,jk), ar.ple.chi2sinit{jk}(q_fiterror) - ...
            ar.ple.chi2s{jk}(q_fiterror), 'r*')        
        plot([0 0]+ar.ple.p(jk), ylim, 'k');
        hold off
        ylabel('\chi^2(init-fit)')
        
        subplot(4,2,5)
        semilogy(ar.ple.ps{jk}(:,jk), sqrt(sum(ar.ple.psinitstep{jk}.^2,2)), '*-')
        hold on
        semilogy(ar.ple.ps{jk}(:,jk), abs(ar.ple.psinitstep{jk}(:,jk)), 'rx-')
        semilogy(xlim, [0 0] + ar.ple.maxstepsize(jk), 'r--')
        semilogy(xlim, [0 0] + ar.ple.minstepsize(jk), 'r--')
        plot([0 0]+ar.ple.p(jk), ylim, 'k');
        hold off
        ylabel('step size')
        
%         subplot(4,2,7)

%         subplot(4,2,[2 4])        
%         plot(ar.ple.ps{jk}(:,jk), ar.ple.betas{jk})
%         hold on
%         plot([0 0]+ar.ple.p(jk), ylim, 'k');
%         hold off
%         ylabel('beta')
% 
%         subplot(4,2,[6 8])
%         plot(ar.ple.ps{jk}(:,jk), ar.ple.alphas{jk})
%         hold on
%         plot([0 0]+ar.ple.p(jk), ylim, 'k');
%         hold off
%         ylabel('alpha')
        
        count = count + 1;
    end
end