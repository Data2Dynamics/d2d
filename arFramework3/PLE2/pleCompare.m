% Compare profile likelihoods of two runs

function pleCompare(ples, labels)

if(nargin==0)
    filenames = fileChooserMulti('./Results', true);
    
    ples = {};
    labels = {};
    for j=1:length(filenames)
        fname = ['./Results/' filenames{j} '/PLE/results.mat'];
        if(exist(fname,'file'))
            tmpple = load(fname);
            ples{end+1} = tmpple.pleGlobals; %#ok<AGROW>
            labels{end+1} = filenames{j}; %#ok<AGROW>
        else
            fprintf('%s does not contains PLE\n', filenames{j});
        end
    end
end

if(isempty(ples))
    return;
end

pLabels = {};
for j=1:length(ples)
    pLabels = union(pLabels, ples{j}.p_labels(ples{j}.q_fit==1));
end

ncols = ceil(length(pLabels)^(0.4))+1;
nrows = ceil(length(pLabels)/ncols);

h = figure(1);
set(h, 'Color', [1 1 1]);

% limits
dchi2 = ples{1}.dchi2_point;
if(ples{1}.plot_simu)
    dchi2 = ples{1}.dchi2;
end

colors = lines(length(ples));
hs = nan(1,length(ples));

for j=1:length(pLabels)
    subplot(nrows,ncols,j);
    
    xvals = [];
    for jple=1:length(ples)
        qj = strmatch(pLabels{j},strvcat(ples{jple}.p_labels),'exact'); %#ok<REMFF1,MATCH3>
        
        if(~isempty(qj))
            % profile
            hs(jple) = plot(ples{jple}.ps{qj}(:,qj), (ples{jple}.chi2s{qj} - ples{jple}.chi2)/dchi2, 'Color', colors(jple,:));
            hold on
            
            % optimum
            plot(ples{jple}.p(qj), (ples{jple}.chi2 - ples{jple}.chi2)/ples{jple}.dchi2, '*', 'Color', colors(jple,:))
            xvals = [xvals; ples{jple}.ps{qj}(:,qj)]; %#ok<AGROW>
        end
    end
    
    xlimtmp2 = (max(xvals)-min(xvals))*0.05;
    if(xlimtmp2>0)
        xlimtmp = [min(xvals)-xlimtmp2 max(xvals)+xlimtmp2];
        xlim(xlimtmp);
    end
    
    % treshold
    plot(xlim, [1 1], 'k--')
    hold off
    
    xlabel(['log_{10}(' strrep(pLabels{j},'_','\_') ')'])
    ylim([-0.1 1.3]);

    if(mod(j-1,ncols)==0)
        ylabel('PL')
    else
        ylabel('')
    end
    
    if(j == length(pLabels))
        legend(hs, strrep(labels,'_','\_'), strrep(labels,'_','\_'));
    end
    set(gca, 'YTickLabel', {})
end


