function [K_org] = motaGetK(percentage,IO)

% otGetK(percentage)

% [K, rowNames] = pwAnalysisOfFits(fits, percentage, gop, threshold1, threshold2)
%
% fits is an array of pwGlobals structs.
%             Default: last applied fit sequence
%
%
% percentage: How many percent of the best fits should be used?
%             Default: 100%
%
% threshold1: Minimum positive correlation for suggested derived parameter with /
%             Default: 0.8
%
% threshold2: Maximum negative correlation for suggested derived parameter with *
%             Default: -0.2
%
% K           (n x p) matrix for the p fitted parameters of the n best fits
%             median is set to 1
%
% K_org       the same as K but with no rescaling
%
%
% rowNames    String array with p names of the fitted parameters
%
%
% Example:    pwAnalysisOfFits([], 20); % Get the best 20% of the last fit sequence
%
%
% Thomas Maiwald, 2006-03-20, 2006-07-26

global pwGlobals2

if~exist('percentage','var')||isempty(percentage)
   percentage=100;
end

if~exist('IO','var')||isempty(IO)
   IO=true;
end

if ~exist('fits','var') || isempty(fits)
    % Get fits of last fitSequence
    if length(pwGlobals2.fitHistory) == 0
        disp('No fits available.')
        return
    end
    lastFitID = pwGlobals2.fitHistory(end).fitID;
    timeOfLastFitSequence = lastFitID(1).timeStart;
    ind = [];
    for i=1:length(pwGlobals2.fitHistory)
        fitID = pwGlobals2.fitHistory(i).fitID;
        if timeOfLastFitSequence == fitID(1).timeStart
            ind = [ind i];
        end
    end
    disp(sprintf('Found %d fits of last fit sequence.', length(ind)));
    fits = pwGlobals2.fitHistory(ind);
end

if length(fits) < 3
    disp(sprintf('Only %d fits available for analysis, but 3 required', length(fits)));
    return
end

% Check whether parsForFitIDs have the same meaning
parsForFitIDs = fits(1).parsForFitIDs(fits(1).indFittedPars);
ok = true;
for i=2:length(fits)
    if ~isequal(parsForFitIDs, fits(i).parsForFitIDs(fits(i).indFittedPars))
        ok = false;
        break;
    end
end
if ok == false
    warndlg('Some fits do not belong to the same fit setting (parsForFitIDs and/or indFittedPars differ).');
    return
end

nFigures = 8;
iFigure  = 1;
for i=1:nFigures
    try
        close(i)
    catch
    end
end
% General_ArrangeFigures(nFigures, pwGlobals2.plotting.menuBar, ...
%     pwGlobals2.plotting.toolBar, 0, pwGlobals2.plotting.distanceToScreenButtom + 150);


filenameTrunk = ['AnalysisOfFits_', datestr(now, 30), '_', fits(1).couples(1).modelID, '_'];

%----------------------------------------------------------------------
% Preparation

% Get matrix fits x parameters
K = zeros(length(fits), length(parsForFitIDs));
for i=1:length(fits)
    K(i,:) = fits(i).parsForFitEnd(fits(i).indFittedPars);
end

K_true = []; % ToDo

% Remove constant parameters, they are
% not fitted and have std 0
stdK = std(K);
if length(find(stdK == 0))>0
    disp(['Some fitted values std = 0!'])
end
K = K(:,find(stdK ~= 0));
K_org = K;
K = K ./ repmat(median(K, 1), size(K, 1), 1);

% Get rownames
rowNames = fits(1).parsForFitIDs(fits(1).indFittedPars);

for i=1:length(rowNames)
    rowNames{i} = regexprep(rowNames{i}, '\s\s+', ' ');
end


%----------------------------------------------------------------------
% Plot chisq values of fits
if IO==true
figure(iFigure)
subplot(1,1,1,'replace')
iFigure = iFigure + 1;
end
chisqValues = [];
for i=1:length(fits)
    chisqValues(i) = fits(i).optimization.results.resNorm;
end
if IO==true
subplot(1,1,1,'replace')
semilogy(chisqValues, 'b.-'); hold on
ylim([0 max(chisqValues) * 1.1])
title(sprintf('Analysis of %d fits (best %g %%)', length(fits), percentage), 'FontSize', 16)
xlabel('fit number', 'FontSize', 14)
ylabel('\chi^2','FontSize', 14)
xlim([0 length(fits)+1])
end
% Use only a certain percentage of fits
if percentage >= 0 && percentage < 100
    nAcceptedFits = fix(percentage * length(chisqValues) / 100);
    [tmp, I] = sort(chisqValues);
    usedFitIDs = I(1:nAcceptedFits);
    if IO==true
    semilogy(usedFitIDs, chisqValues(usedFitIDs), 'ro');    
    end
    K = K(usedFitIDs, :);
    K_org = K_org(usedFitIDs, :);
    K = K ./ repmat(median(K, 1), size(K, 1), 1);    
end
