% arCompareWithBiobaseSimulation(file)
% 
%  Reads the download of Biobase -> Actions -> BioModels Online Simulation
%  and compares with ar.model(1).condition(1).xFineSimu
% 
%  Be sure that xFineSimu has been simulated for the current parameters,
%  e.g. by calling
%       arQplot('x');
%       arPlot
% 
% Example:
%   arQplot('x');
%   arPlot
%   file = 'SIMU1447940350922.dat'
%   arCompareWithBiobaseSimulation(file)

function arCompareWithBiobaseSimulation(file)

global ar

m = 1;
c = 1;

dat = importdata(file);

dat.x = strsplit(dat.textdata{1},' ');
dat.indcol = find(sum(dat.data~=0,1)>0);

t = dat.data(:,1);

subx = ceil(sqrt(length(dat.indcol)-1));
suby = ceil((length(dat.indcol)-1)/subx);

for i=2:length(dat.x)
    if(length(dat.x{i})==1)
        dat.x{i} = [dat.x{i},'_state'];
    end
end

if(length(dat.indcol)<16)
    fs = 10;
elseif(length(dat.indcol)<25)
    fs = 8;
elseif(length(dat.indcol)<36)
    fs = 7;
else
    fs = 6;
end


%%
close all
dolegend = 1;
for i=length(dat.indcol):-1:2
    subplot(subx,suby,i-1)
    set(gca,'FontSize',fs)
    plot(t,dat.data(:,dat.indcol(i)),'k');
    hold on
    ind = strmatch(dat.x{dat.indcol(i)},ar.model(m).x,'exact');
    
    if isempty(ind)
        if isempty(strmatch(dat.x{dat.indcol(i)},ar.pLabel,'exact'));
            fprintf('%s from BIOMODELs simulation neither found as dynamic state nor as parameter.\n',dat.x{dat.indcol(i)})
        end
    elseif(length(ind)>1)
        warning(sprintf('%s from BIOMODELs simulation multiple times.\n',dat.x{dat.indcol(i)}))
    else
        plot(t,interp1(ar.model(m).condition(c).tFine,ar.model(m).condition(c).xFineSimu(:,ind),t),'r--')
        if dolegend ==1
            set(legend('Biobase','d2d'),'FontSize',fs);
            dolegend = 0; % only once
        end
    end    
    xlim([0,100])
    title(strrep(dat.x{dat.indcol(i)},'_','\_'),'FontSize',fs)
end


saveas(gcf,'arCompareWithBiobaseSimulation');


