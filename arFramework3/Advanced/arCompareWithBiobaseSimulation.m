% pass = arCompareWithBiobaseSimulation(file, [silent])
% 
% Reads the download of Biobase -> Actions -> BioModels Online Simulation
% and compares with ar.model(1).condition(1).xFineSimu
%
%   file   - file from database simulation, eg 'SIMU*.dat'
%   silent - boolean for plotting (ar and biomodels in same plot) [false]
%
%   pass   - boolean, if relative difference <rtol, pass = 1 
% 
% Be sure that xFineSimu has been simulated for the current parameters,
% e.g. by calling
%      arQplot('x');
%      arPlot
% 
% Example:
%   arQplot('x');
%   arPlot
%   file = 'SIMU1447940350922.dat'
%   arCompareWithBiobaseSimulation(file)

function pass = arCompareWithBiobaseSimulation(file, silent)

global ar

m = 1;
c = 1;

% Tolerances for a pass
rtol = 1e-3;
atol = 1e-4;

if (nargin<2)
    silent = false;
end

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
    biomodelsSim = dat.data(:,dat.indcol(i));
    if ( ~silent )
        subplot(subx,suby,i-1)
        set(gca,'FontSize',fs)
        plot(t, biomodelsSim,'k');
        hold on
    end
    
    ind = strmatch(dat.x{dat.indcol(i)},ar.model(m).x, 'exact'); %#ok
    if isempty(ind)
        if isempty(strmatch(dat.x{dat.indcol(i)},ar.pLabel, 'exact')); %#ok
            arFprintf(2, '%s from BIOMODELs simulation neither found as dynamic state nor as parameter.\n',dat.x{dat.indcol(i)})
        end
    elseif(length(ind)>1)
        warning('%s from BIOMODELs simulation multiple times.\n',dat.x{dat.indcol(i)})
    else
        d2dSim           = interp1(ar.model(m).condition(c).tFine,ar.model(m).condition(c).xFineSimu(:,ind),t);
        d2dFilt          = bsxfun(@max, d2dSim, atol);
        biomodelsFilt    = bsxfun(@max, biomodelsSim, atol);
        maxDifference(i) = max( ( (d2dFilt - biomodelsFilt) ./ (biomodelsFilt) ).^2 );

        if ( ~silent )
            plot(t, d2dSim,'r--')
            if dolegend ==1
                set(legend('Biobase','d2d'),'FontSize',fs);
                dolegend = 0; % only once
            end
        end
    end
    
    if ( ~silent )
        xlim([0,100]);
        title(strrep(dat.x{dat.indcol(i)},'_','\_'),'FontSize',fs);
    end
end

% Fit acceptable?
if ( max( maxDifference ) > rtol )
    pass = 0;
else
    pass = 1;
end

if ( ~silent )
    saveas(gcf,'arCompareWithBiobaseSimulation');
end


