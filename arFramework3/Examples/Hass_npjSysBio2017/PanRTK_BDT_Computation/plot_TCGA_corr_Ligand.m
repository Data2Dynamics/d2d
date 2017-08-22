%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   plot_TCGA_corr_Ligand                   %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  script to plot TCGA patient-data correlation between ligand expression
%  and predicted response for breast/colorectal cancer
% Plotting is done via Violin plots via a script obtained from:
% https://de.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions--distributionplot-m-


function plot_TCGA_corr_Ligand()
global TCGA_Matrices;
fields = fieldnames(TCGA_Matrices);

figure
nr_TCGA=sum(~cellfun(@isempty,strfind(fields,'TCGA_')));   
 ligands = {'EGF','HRG','HGF'};
 ifig=0;
for i=1:length(fields)
    if(isempty(strfind(fields{i},'TCGA_')))
        continue;
    end
    ifig = ifig+1;
    lig = TCGA_Matrices.(fields{i});
    for j=1:length(ligands)
        outcome = TCGA_Matrices.([strrep(fields{i},'TCGA_','') '_' ligands{j}]);
%         outcome_prob = TCGA_Matrices.([strrep(fields{i},'TCGA_','') '_' ligands{j} '_Prob']);        
        neg=lig(5+j,strcmp(outcome,'0'));
        pos=lig(5+j,strcmp(outcome,'1'));       
        
        [h, p_val] = ttest2(neg,pos,'Alpha',0.05);
%        p_val
        if(length(neg)<length(pos))
           neg = [neg NaN(1,length(pos)-length(neg))]; 
        else
           pos = [pos NaN(1,length(neg)-length(pos))]; 
        end

        mat = [neg' pos'];
        
       
        subplot(5,3,(ifig-1)*length(ligands)+j)
%         boxplot(log(mat+1),'Labels',{'non-responder','responder'},'PlotStyle','compact','notch','on','colorgroup',{'non-responder','responder'},'colors',[0         0    0.5625;0.5000         0         0])
        distributionPlot(log2(mat+1),'showMM',6,'color',{'dg','dr'},'xNames',{'',''})%[0         0    0.5625;0.5000         0         0]) ,'xNames',{'non-responder','responder'}
        set(gcf,'Color','w')
        ylabel([ligands{j} ' log2-expression'])
        %title({strrep(fields{i},'TCGA_',''); ['Response for ' ligands{j}]})
        title(['Response for ' ligands{j}])
        h1 = get(gca,'Children');
        min_tmp = 0;
        max_tmp = 0;
        for jh=1:length(h1)
            if(strcmp(h1(jh).Type,'patch'))
                if(min(h1(jh).YData(1,:))<min_tmp)
                    min_tmp = min(h1(jh).YData(1,:));
                end
                if(max(h1(jh).YData(1,:))>max_tmp)
                    max_tmp = max(h1(jh).YData(1,:));                   
                end
            end
        end       
        set(gca,'YLim',[[min_tmp max_tmp] + [0 2]])
        YLim_tmp=get(gca,'YLim');   
        text2 ='-fold';
        if(h)
            if((p_val>0.01 || (ifig==1 && j==1)))
                text2 = [text2 '(*)'];
            elseif(ifig~=2 && j~=1)
                text2 = [text2 '(**)'];
            end
        end
        if(ifig==2 && j==1)
            text2 = [text2 '(**)'];
        end
        text(.7,YLim_tmp(2)-0.5,[num2str(round(2^(nanmean(log2(pos+1))-nanmean(log2(neg+1))),1)) text2])
        
        set(gca,'XTickLabelRotation',35);
        save_folder='./TCGA_Predictions/';
        set(gcf,'PaperSize',[20 22]);

        saveas(gcf,[save_folder 'TCGA_ligand_correlation.fig'])
        saveas(gcf,[save_folder 'TCGA_ligand_correlation'],'svg') 
        
    end
end