%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   dose_plots                              %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Do plotting of signaling model for every ligand, cell line and
%  co-stimulations separately
%
%
function dose_plots(save, which)
global ar
if(~exist('save','var'))
    save = false;
end

if(~exist('which','var') || isempty(which))
    which=1:length(ar.model(1).plot);
end
ar.model(1).qPlotXs(:) = 0;
arSimu(false,true,true)
for ip = which
    nd = length(ar.model(1).plot(ip).dLink);
    ar.model(1).qPlotYs(:) = 0;
    ar.model(1).qPlotYs(ip) = 1;
    conditions = NaN(length(ar.model(1).data(1).condition),nd);
    for j=1:nd
        for i=1:length(ar.model(1).data(ar.model(1).plot(ip).dLink(j)).condition)
            conditions(i,j)=str2double(ar.model(1).data(ar.model(1).plot(ip).dLink(j)).condition(i).value);
        end   
    end

    axis_high = [];
    axis_low = [];
    
    for i = 1:nd
        ac_data = ar.model(1).plot(ip).dLink(i);
        qUnlog = ar.model(1).data(ac_data).logfitting & ...
                        ~ar.model(1).data(ac_data).logplotting;
        if(sum(qUnlog)>size(ar.model(1).data(ac_data).yExp,2)/2)           
            axis_high = max([axis_high; 10.^(ar.model(1).data(ac_data).yFineSimu+ar.model(1).data(ac_data).ystdFineSimu); 10.^ar.model(1).data(ac_data).yExp]);
            axis_low = min([axis_low; 10.^(ar.model(1).data(ac_data).yFineSimu-ar.model(1).data(ac_data).ystdFineSimu); 10.^ar.model(1).data(ac_data).yExp]);
        else
            axis_high = max([axis_high; ar.model(1).data(ac_data).yFineSimu+ar.model(1).data(ac_data).ystdFineSimu; ar.model(1).data(ac_data).yExp]);
            axis_low = min([axis_low; ar.model(1).data(ac_data).yFineSimu-ar.model(1).data(ac_data).ystdFineSimu; ar.model(1).data(ac_data).yExp]);
        end
    end

    EGFs = [];
    IGFs = [];
    HRGs = [];
    HGFs = [];
    BTCs = [];
    cos = [];
    conds = {'EGF_level','HRG_level','HGF_level','IGF1_level','BTC_level'};
    conds_id = NaN(1,length(conds));
    for j=1:length(ar.model(1).data(ar.model(1).plot(ip).dLink(1)).condition)
        idC = ismember(conds,ar.model(1).data(ar.model(1).plot(ip).dLink(1)).condition(j).parameter);
        conds_id(find(idC==1)) = j;
    end
    for j=1:nd
        condID = conds_id(~isnan(conds_id));
        ac_data = ar.model(1).plot(ip).dLink(j);
        getZero=conditions(condID,j)~=0;
        if(sum(getZero)==0)
            EGFs = [EGFs j]; IGFs = [IGFs j]; HRGs = [HRGs j]; HGFs = [HGFs j]; BTCs = [BTCs j]; cos = [cos j];
        elseif(sum(getZero)==1)
            hae = condID(find(conditions(condID,j)>0));
            if(strcmp(conds{find(conds_id==hae)},'EGF_level'))
                EGFs = [EGFs j];
            elseif(strcmp(conds{find(conds_id==hae)},'HRG_level'))
                HRGs = [HRGs j]; 
            elseif(strcmp(conds{find(conds_id==hae)},'IGF1_level'))
                IGFs = [IGFs j];
            elseif(strcmp(conds{find(conds_id==hae)},'HGF_level'))
                HGFs = [HGFs j];
            elseif(strcmp(conds{find(conds_id==hae)},'BTC_level'))
                BTCs = [BTCs j];
            end
        end
        if((sum(getZero)==2 && conditions(conds_id(4),j)==0) || (sum(getZero)==1 && sum(conditions(conds_id(~isnan(conds_id)),j))==2.5))
            cos = [cos j];
        end
    end
    
    down_plot = {'pERK_au','pAKT_au','pS6_au','pMEK_au'};
    
    plot_range = [];
    for i=1:length(down_plot)
         plot_range = [plot_range strmatch(down_plot{i},ar.model(1).yNames)+1];
    end
    
    ar.model(1).plot(ip).cond_bkp = ar.model(1).plot(ip).condition;
    ar.model(1).plot(ip).dLink_bkp = ar.model(1).plot(ip).dLink;
        
    if(length(EGFs)>1)
        ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp(EGFs);
        ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp(EGFs);
        arPlot(0,0,0,0)
        curName = get(gcf,'Name');
    
        newName=['1' curName];
        set(gcf,'Name',newName');
        h2 = get(gcf,'Children');
        for i=1:length(h2)
            if(strcmp(h2(i).Type,'axes'))
                title = strrep(h2(i).Title.String,'\_au','');
                if(isempty(title))
                    continue
                end
                yindex = find(strcmp(strrep(ar.model(1).data(1).yNames,'_au',''),title{1}));
                if(strcmp(title,'pIGF1R'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pIGF1R','pIGF-1R');
                elseif(strcmp(title,'pErbB2'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pErbB2','pHER2');
                end
                if(~isempty(yindex))
                    if(axis_low(yindex)>0)
                        set(h2(i),'YLim',[axis_low(yindex)*.9 axis_high(yindex)*1.1]);
                    else
                        set(h2(i),'YLim',[axis_low(yindex)*1.1 axis_high(yindex)]);
                    end
                end
            end
        end
        plot_range_tmp = [];
        RTK_plot = {'pEGFR_au','pErbB2_au','pErbB3_au'};

        for i=1:length(RTK_plot)
             plot_range_tmp = [plot_range_tmp strmatch(RTK_plot{i},ar.model(1).yNames)+1];
        end
        plot_range_tmp = [plot_range_tmp plot_range];
        if(save)
            saveName = sprintf('EGF_%i',ip);
            hhRearrangeSubplots(plot_range_tmp,1,1,saveName,gcf)
            %arSaveFigure(gcf,[strrep(curName,'Y: ',''),'_EGF'],'/Figures/Y')
        else
            %hhRearrangeSubplots([10 11 12 plot_range],0,1,[],gcf)
        end
    end

    if(length(HRGs)>1)
        %set(gcf,'NextPlot','new');
        ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp(HRGs);
        ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp(HRGs);
        arPlot(0,0,0,0)
        newName=['2' curName];
        set(gcf,'Name',newName');
        h2 = get(gcf,'Children');
        for i=1:length(h2)
            if(strcmp(h2(i).Type,'axes'))
                title = strrep(h2(i).Title.String,'\','');
                if(isempty(title))
                    continue
                end
                yindex = find(strcmp(strrep(ar.model(1).data(1).yNames,'_au',''),title{1}));
                if(strcmp(title,'pIGF1R'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pIGF1R','pIGF-1R');
                elseif(strcmp(title,'pErbB2'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pErbB2','pHER2');
                end
                if(~isempty(yindex))
                    if(axis_low(yindex)>0)
                        set(h2(i),'YLim',[axis_low(yindex)*.9 axis_high(yindex)*1.1]);
                    else
                        set(h2(i),'YLim',[axis_low(yindex)*1.1 axis_high(yindex)]);
                    end
                end
            end
        end
        plot_range_tmp = [];
        RTK_plot = {'pEGFR_au','pErbB2_au','pErbB3_au'};

        for i=1:length(RTK_plot)
             plot_range_tmp = [plot_range_tmp strmatch(RTK_plot{i},ar.model(1).yNames)+1];
        end
        plot_range_tmp = [plot_range_tmp plot_range];
        if(save)
            saveName = sprintf('HRG_%i',ip);
            %arSaveFigure(gcf,[strrep(curName,'Y: ',''),'_IGF'],'/Figures/Y')
            hhRearrangeSubplots(plot_range_tmp,1,1,saveName,gcf)
        else
            %hhRearrangeSubplots([11 12 plot_range],0,1,[],gcf)
        end
    end
    if(length(IGFs)>1)
        ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp(IGFs);
        ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp(IGFs);
        arPlot(0,0,0,0)
        newName=['3' curName];
        set(gcf,'Name',newName');
        h2 = get(gcf,'Children');
        for i=1:length(h2)
            if(strcmp(h2(i).Type,'axes'))
                title = strrep(h2(i).Title.String,'\','');
                if(isempty(title))
                    continue
                end
                yindex = find(strcmp(strrep(ar.model(1).data(1).yNames,'_au',''),title{1}));
                if(strcmp(title,'pIGF1R'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pIGF1R','pIGF-1R');
                elseif(strcmp(title,'pErbB2'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pErbB2','pHER2');
                end
                if(~isempty(yindex))
                    if(axis_low(yindex)>0)
                        set(h2(i),'YLim',[axis_low(yindex)*.9 axis_high(yindex)*1.1]);
                    else
                        set(h2(i),'YLim',[axis_low(yindex)*1.1 axis_high(yindex)]);
                    end
                end
            end
        end
        plot_range_tmp = [];
        RTK_plot = {'pIGF1R_au','pEGFR_au','pErbB2_au'};

        for i=1:length(RTK_plot)
             plot_range_tmp = [plot_range_tmp strmatch(RTK_plot{i},ar.model(1).yNames)+1];
        end
        plot_range_tmp = [plot_range_tmp plot_range];
        if(save)
            saveName = sprintf('IGF1_%i',ip);
            hhRearrangeSubplots(plot_range_tmp,1,1,saveName,gcf)
            %arSaveFigure(gcf,[strrep(curName,'Y: ',''),'_HRG'],'/Figures/Y')
        else
            %hhRearrangeSubplots([11 13 plot_range],0,1,[],gcf)
        end
    end
    if(length(HGFs)>1)
        ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp(HGFs);
        ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp(HGFs);
        arPlot(0,0,0,0)
        newName=['4' curName];
        set(gcf,'Name',newName');
        h2 = get(gcf,'Children');
        for i=1:length(h2)
            if(strcmp(h2(i).Type,'axes'))
                title = strrep(h2(i).Title.String,'\','');
                if(isempty(title))
                    continue
                end
                yindex = find(strcmp(strrep(ar.model(1).data(1).yNames,'_au',''),title{1}));
                if(strcmp(title,'pIGF1R'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pIGF1R','pIGF-1R');
                elseif(strcmp(title,'pErbB2'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pErbB2','pHER2');
                end
                if(~isempty(yindex))
                    if(axis_low(yindex)>0)
                        set(h2(i),'YLim',[axis_low(yindex)*.9 axis_high(yindex)*1.1]);
                    else
                        set(h2(i),'YLim',[axis_low(yindex)*1.1 axis_high(yindex)]);
                    end
                end
            end
        end
        plot_range_tmp = [];
        RTK_plot = {'pEGFR_au','pErbB3_au','pMet_au'};

        for i=1:length(RTK_plot)
             plot_range_tmp = [plot_range_tmp strmatch(RTK_plot{i},ar.model(1).yNames)+1];
        end
        plot_range_tmp = [plot_range_tmp plot_range];
        if(save)
            saveName = sprintf('HGF_%i',ip);
            hhRearrangeSubplots(plot_range_tmp,1,1,saveName,gcf)
            %arSaveFigure(gcf,[strrep(curName,'Y: ',''),'_HRG'],'/Figures/Y')
        else
            %hhRearrangeSubplots([11 13 plot_range],0,1,[],gcf)
        end
    end
    if(length(BTCs)>1)
        ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp(BTCs);
        ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp(BTCs);
        arPlot(0,0,0,0)
        newName=['5' curName];
        set(gcf,'Name',newName');
        h2 = get(gcf,'Children');
        for i=1:length(h2)
            if(strcmp(h2(i).Type,'axes'))
                title = strrep(h2(i).Title.String,'\','');
                if(isempty(title))
                    continue
                end
                yindex = find(strcmp(strrep(ar.model(1).data(1).yNames,'_au',''),title{1}));
                if(strcmp(title,'pIGF1R'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pIGF1R','pIGF-1R');
                elseif(strcmp(title,'pErbB2'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pErbB2','pHER2');
                end
                if(~isempty(yindex))
                    if(axis_low(yindex)>0)
                        set(h2(i),'YLim',[axis_low(yindex)*.9 axis_high(yindex)*1.1]);
                    else
                        set(h2(i),'YLim',[axis_low(yindex)*1.1 axis_high(yindex)]);
                    end
                end
            end
        end
        plot_range_tmp = [];
        RTK_plot = {'pEGFR_au','pErbB2_au','pErbB3_au'};

        for i=1:length(RTK_plot)
             plot_range_tmp = [plot_range_tmp strmatch(RTK_plot{i},ar.model(1).yNames)+1];
        end
        plot_range_tmp = [plot_range_tmp plot_range];
        if(save)
            saveName = sprintf('BTC_%i',ip);
            hhRearrangeSubplots(plot_range_tmp,1,1,saveName,gcf)
            %arSaveFigure(gcf,[strrep(curName,'Y: ',''),'_HRG'],'/Figures/Y')
        else
            %hhRearrangeSubplots([11 13 plot_range],0,1,[],gcf)
        end
    end
    if(length(cos)>3)
        ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp(cos);
        ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp(cos);
        if(ip==10)
           ar.model(1).plot(ip).condition =  ar.model(1).plot(ip).condition([1 5 12 14]);
           ar.model(1).plot(ip).dLink =  ar.model(1).plot(ip).dLink([1 5 12 14]);
        end
        arPlot(0,0,0,0)
        newName=['5' curName];
        set(gcf,'Name',newName');
        h2 = get(gcf,'Children');
        for i=1:length(h2)
            if(strcmp(h2(i).Type,'axes'))
                title = strrep(h2(i).Title.String,'\','');
                if(isempty(title))
                    continue
                end
                yindex = find(strcmp(strrep(ar.model(1).data(1).yNames,'_au',''),title{1}));
                if(strcmp(title,'pIGF1R'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pIGF1R','pIGF-1R');
                elseif(strcmp(title,'pErbB2'))
                    h2(i).Title.String = strrep(h2(i).Title.String,'pErbB2','pHER2');
                end
                if(~isempty(yindex))
                    if(axis_low(yindex) > 0)
                        set(h2(i),'YLim',[axis_low(yindex)*.9 axis_high(yindex)*1.1]);
                    else
                        set(h2(i),'YLim',[axis_low(yindex)*1.1 axis_high(yindex)]);
                    end
                end
            end
        end
        plot_range_tmp = [];
        RTK_plot = {'pEGFR_au','pErbB2_au','pErbB3_au','pIGF1R_au','pMet_au'};

        for i=1:length(RTK_plot)
             plot_range_tmp = [plot_range_tmp strmatch(RTK_plot{i},ar.model(1).yNames)+1];
        end
        plot_range_tmp = [plot_range_tmp plot_range];
        if(save)
            saveName = sprintf('coStim_%i',ip);
            hhRearrangeSubplots(plot_range_tmp,1,1,saveName,gcf)
            %arSaveFigure(gcf,[strrep(curName,'Y: ',''),'_cos'],'/Figures/Y')
        else
           %hhRearrangeSubplots([10:13 plot_range],0,1,[],gcf)
        end
    end
    ar.model(1).plot(ip).condition = ar.model(1).plot(ip).cond_bkp;
    ar.model(1).plot(ip).dLink = ar.model(1).plot(ip).dLink_bkp;
    if(save)
        close all
    end

end