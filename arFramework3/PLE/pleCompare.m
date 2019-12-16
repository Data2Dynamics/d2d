% pleCompare([ples], [labels], [savetofile])
% 
% Compare profile likelihoods of two runs
%
%    ples              cell array with ple structs to compare, 
%                      if is empty fileChooser is opened
%    labels            cell array with labels for the ples,
%                      if is empty fileChooser is opened
%    savetofile  [0]   boolean, specifies if plot with results is saved to
%                      [ples{end}.savePath '/ple_compare']. 

function pleCompare(ples, labels, savetofile)

if(nargin==0) || (isempty(ples) || isempty(labels))
    filenames = fileChooserMulti('./Results', true);
    
    ples = {};
    labels = {};
    for j=1:length(filenames)
%         fname = ['./Results/' filenames{j} '/PLE/results.mat']
        fname = ['./Results/' filenames{j} '/workspace.mat']
        if(exist(fname,'file'))
            tmpple = load(fname);
            ples{end+1} = tmpple.ar.ple; %#ok<AGROW>
            labels{end+1} = filenames{j}; %#ok<AGROW>
        else
            fprintf('%s does not contains PLE\n', filenames{j});
        end
        ples{end}.chi2 = min(cell2mat(cellfun(@min,ples{end}.chi2s, 'UniformOutput', false)));
    end
end

drin = [];
for i=1:length(ples)
    if(~isempty(ples{i}.ps))
        drin = [drin,i];
    end
end
if(isempty(drin))
    warning('pleCompare.m: No calculated ples provided!')
    return
else
    ples = ples(drin);
    labels = labels(drin);
end


pLabels = {};
for j=1:length(ples)
%     qq = ples{j}.q_fit==1;
%     qq(1:length(ples{j}.ps)) = qq(1:length(ples{j}.ps)) & ~cellfun(@isempty, ples{j}.ps);
    qq =  ~cellfun(@isempty, ples{j}.ps);
    qq((length(ples{j}.ps)+1):end) = false;
    pLabels = union(pLabels, ples{j}.p_labels(qq)); %R2013a compatible
end

ncols = ceil((length(pLabels)+1)^(0.4))+1;
nrows = ceil((length(pLabels)+1)/ncols);

h = figure(1);
set(h, 'Color', [1 1 1]);

% limits
if(isfield(ples{1}, 'dchi2_point'))
    dchi2 = ples{1}.dchi2_point;
    if(ples{1}.plot_simu)
        dchi2 = ples{1}.dchi2;
    end
elseif(isfield(ples{1}, 'alpha') && isfield(ples{1}, 'ndof'))
    dchi2 = arChi2inv(1-ples{1}.alpha, ples{1}.ndof);
else
    error('no information on dchi2, alpha or ndof');
end


colors = lines(length(ples));
hs = nan(1,length(ples));

for j=1:length(pLabels)
    subplot(nrows,ncols,j);
    
    xvals = [];
    for jple=1:length(ples)
        qj = strmatch(pLabels{j},strvcat(ples{jple}.p_labels),'exact'); %#ok<REMFF1,MATCH3>
       %  
        if(~isempty(qj) && size(ples{jple}.ps,2)>=qj && ~isempty(ples{jple}.ps{qj}))
            % profile
            lw = 1;
            marker = '-.';
            hs(jple) = plot(ples{jple}.ps{qj}(:,qj), (ples{jple}.chi2s{qj} - ples{jple}.chi2)/dchi2, marker,'Color', colors(jple,:), ...
                'LineWidth', lw);
            hold on
            
            % optimum
            plot(ples{jple}.p(qj), (ples{jple}.chi2 - ples{jple}.chi2)/dchi2, '*', 'Color', colors(jple,:))
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
    set(gca, 'YTickLabel', {})
    
    if(j == length(pLabels))
        subplot(nrows,ncols,j+1)
        mypos = get(gca,'Position');
        myleg = legend(hs, strrep(labels,'_','\_'), strrep(labels,'_','\_'));
        legpos = get(myleg,'Position');
        set(myleg,'Position',[mypos(1) mypos(2)+mypos(4)-legpos(4) legpos(3:4)])
        set(gca,'Visible','off')
    end
end

if ~exist('savetofile','var')
    savetofile = false;
end

% save
if(savetofile && exist(ples{end}.savePath, 'dir'))
    saveas(gcf, [ples{end}.savePath '/ple_compare'], 'fig')
    print('-depsc2', [ples{end}.savePath '/ple_compare']);
    if(ispc)
        print('-dpdf', [ples{end}.savePath '/ple_compare']);
    else
        system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' [ples{end}.savePath '/ple_compare'] '.eps '  [ples{end}.savePath '/ple_compare'] '.pdf']);
    end
end


