%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   plot_BDT                                %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting script for ROC curves
% struct = the input struct in the global 'bdt_figures', i.e. the save name
% you put in BDT_bootstrap

% folder = output folder for the figures, if blank, BDT_plots is used

% which_lig = which ligand to be plotted, same as for BDT_bootstrap

function plot_BDT(figure_struct, folder, which_lig)
global bdt_figures
if(~exist('figure_struct','var') || isempty(figure_struct))
    error('please provide the Classifier');
elseif(ischar(figure_struct))
    struct_tmp = {figure_struct};
else
    struct_tmp = figure_struct;
end

if(~isfield(bdt_figures,figure_struct))
    error('cannot find the field in bdt_figures struct');
end

if(~exist('folder','var') || isempty(folder))
    folder='BDT_plots';
end

% check if dir exists
if(~isempty(folder) && ~exist(folder, 'dir'))
    mkdir(folder)
end

%find right struct for single MM/Ligand tests
which_drug = 'CONTROL';

if(~exist('which_lig','var') || isempty(which_lig) || (isinteger(which_lig) && which_lig==1))
    which_lig='ALL';
elseif((isinteger(which_lig) && which_lig==0))
    which_lig = 'CONTROL';
else
    which_lig = upper(which_lig );
end

if(~isempty(which_lig) && isempty(which_drug))
    which_drug='ALL';
elseif(~isempty(which_drug) && isempty(which_lig))
    which_lig='ALL';
end

if(~ischar(which_drug) || ~ischar(which_lig))
    error('please specify correct drug/ligand strings');
end

plot_id=[];

%for all inputs, do ROC curve plot and save it
for i=1:length(struct_tmp)   
    if(~isempty(which_drug))
       for j=1:length(bdt_figures.(struct_tmp{i})(1).testing)
            if(strcmp(which_drug,bdt_figures.(struct_tmp{i})(1).testing(j).which_drugs) && strcmp(which_lig,bdt_figures.(struct_tmp{i})(1).testing(j).which_ligs))
               plot_id = j; 
            end
       end
    end
    if(isempty(plot_id))
        error('couldnt find the specified combination in the struct');
    end
    save_name1 = ['ROC_' bdt_figures.(struct_tmp{i})(1).name];
    save_name1 = ['./' folder '/' save_name1];
    if(~isempty(plot_id))
        save_name1= [save_name1 '_drug_' strrep(which_drug,'/',' and ') '_lig_' strrep(which_lig,'/',' and ') ];
        
    end
    titles={'Significant shrinkage','No sign. proliferation','Significant proliferation'};
    save_names={'shrinkage','neutral','proliferation'};
    if(strfind(bdt_figures.(struct_tmp{i})(1).name,'RNA'))
        title_name = 'RNA, ';
    else
        title_name = 'model, ';
    end
    if(~isempty(plot_id))
        title_name = [title_name 'drug ' strrep(which_drug,'_','\_') ', ligand ' which_lig];
    end
    
    uXs = [];
    uClass =[];
    uYs = [];
    %get unique classes/length of sets
    for j=1:length(bdt_figures.(struct_tmp{i}))
        if(size(bdt_figures.(struct_tmp{i})(j).testing(plot_id).xVal,2)>length(uXs))
            uXs = size(bdt_figures.(struct_tmp{i})(j).testing(plot_id).xVal,2);
            uClass =  bdt_figures.(struct_tmp{i})(j).testing(plot_id).bdt_class;
            uYs = bdt_figures.(struct_tmp{i})(j).testing(plot_id).Y_unique;
        end
    end
    %for all outcomes there is an individual roc curve    
    for j=1:uXs
        
        if(length(uClass)<j || length(uYs)<j)
            if(isempty(uClass))
                error('Apparently there is no data in the condition you wanted to plot, do a meaningful testing with re-use in the BDT script');
            end
           continue; 
        end
        struct_unique=struct();
        auc_norm = [];
        auc_rnd = [];
        
        x_all_norm = [];
        y_all_norm = [];
        x_all_rnd = [];
        y_all_rnd = [];
        
        class = uClass(j);
        title_tmp = titles{class+2};
        save_name = [save_name1 save_names{class+2}];
        %Gather unique X-Y values of all trees, change here if you need a
        %different X,Y from struct
        for k=1:length(bdt_figures.(struct_tmp{i}))
            which_class = find(bdt_figures.(struct_tmp{i})(k).testing(plot_id).bdt_class==class);
             if(isempty(which_class))
                continue;
            end
            x_tmp = unique(bdt_figures.(struct_tmp{i})(k).testing(plot_id).xVal(:,which_class));
           
            x_tmp(isnan(x_tmp))=[];
            y_tmp = [];
            if(isempty(x_tmp))
               continue; 
            end
            for l=1:length(x_tmp)
                %gather unique X and Y values for each ROC curve
                where = find(bdt_figures.(struct_tmp{i})(k).testing(plot_id).xVal(:,which_class)==x_tmp(l));
                y_tmp(l) = bdt_figures.(struct_tmp{i})(k).testing(plot_id).yVal(where(end),which_class);
                
            end           
            struct_unique(k).yVal = y_tmp';
            struct_unique(k).xVal = x_tmp;
            x_all_norm = [x_all_norm x_tmp'];
            
            %put together all area under curve measures of ROC curves
            auc_norm = [auc_norm bdt_figures.(struct_tmp{i})(k).testing(plot_id).auc(which_class)];
            
            if(isfield(bdt_figures.(struct_tmp{i})(k),'tree_rnd') && ~isempty(bdt_figures.(struct_tmp{i})(k).tree_rnd))
                %do the same for randomized outcome sample
                which_rnd_class = find(bdt_figures.(struct_tmp{i})(k).testing(plot_id).bdt_rnd_class==class);            
                if(isempty(which_rnd_class))
                    continue;
                end
                x_tmp = unique(bdt_figures.(struct_tmp{i})(k).testing(plot_id).xVal_rnd(:,which_rnd_class));
                x_tmp(isnan(x_tmp))=[];
                y_tmp=[];

                if(isempty(x_tmp))
                   continue; 
                end
                for l=1:length(x_tmp)
                    where = find(bdt_figures.(struct_tmp{i})(k).testing(plot_id).xVal_rnd(:,which_rnd_class)==x_tmp(l));
                    y_tmp(l) = bdt_figures.(struct_tmp{i})(k).testing(plot_id).yVal_rnd(where(end),which_rnd_class);               
                end
                struct_unique(k).yVal_rnd = y_tmp';
                struct_unique(k).xVal_rnd = x_tmp;
                x_all_rnd = [x_all_rnd x_tmp'];
                auc_rnd = [auc_rnd bdt_figures.(struct_tmp{i})(k).testing(plot_id).auc_rnd(which_rnd_class)];
            end
            
        end      
        if(isempty(x_all_norm) || isempty(x_all_rnd))
           sprintf('NOT ENOUGH EVENTS IN CATEGORY %s , WILL CONTINUE WITH NEXT \n',title_tmp); 
        end
        %get unique x's, interpolate all ROC curves to these
%         x_norm_u = unique(sort(x_all_norm));
%         x_rnd_u = unique(sort(x_all_rnd));
        x_norm_u = 0:0.01:1;
        x_rnd_u = 0:0.01:1;
        y_tmp = [];
        y_rnd_tmp = [];
        for k=1:length(struct_unique)            
            if(length(struct_unique(k).xVal)>0)
                bdt_figures.(struct_tmp{i})(k).testing(plot_id).interp(j).yVal = interp1(struct_unique(k).xVal,struct_unique(k).yVal,x_norm_u,'nearest','extrap');  %'nearest'         
                if(isnan(y_tmp))
                    y_tmp = [];
                end
                y_tmp(:,k) = bdt_figures.(struct_tmp{i})(k).testing(plot_id).interp(j).yVal;
            else
                y_tmp(:,k) = NaN;
            end
            if(isfield(struct_unique(k),'xVal_rnd') && length(struct_unique(k).xVal_rnd)>0)
                if(isnan(y_rnd_tmp))
                    y_rnd_tmp = [];
                end
                bdt_figures.(struct_tmp{i})(k).testing(plot_id).interp(j).yVal_rnd = interp1(struct_unique(k).xVal_rnd,struct_unique(k).yVal_rnd,x_rnd_u,'nearest','extrap');
                y_rnd_tmp(:,k) = bdt_figures.(struct_tmp{i})(k).testing(plot_id).interp(j).yVal_rnd;
            else
                y_rnd_tmp(:,k) = NaN;
            end            
        end
        
        %save all mean/std/auc values in struct
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).mean = nanmean(y_tmp,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_mean = nanmean(y_rnd_tmp,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).std = nanstd(y_tmp,1,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_std = nanstd(y_rnd_tmp,1,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).x_unique = x_norm_u;
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).x_rnd_unique = x_rnd_u;
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).auc_norm = auc_norm;
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).auc_rnd = auc_rnd;
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).auc_mean = round(100*mean(auc_norm),2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).auc_rnd_mean = round(100*mean(auc_rnd),2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).auc_std = round(100*std(auc_norm),2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).auc_rnd_std = round(100*std(auc_rnd),2);
        
        %do patches
        tmpx_norm = [x_norm_u'; flipud(x_norm_u')];
        tmpy_norm = [bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).mean+bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).std; flipud(bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).mean-bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).std)];
        tmpx_rnd = [x_rnd_u'; flipud(x_rnd_u')];
        tmpy_rnd = [bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_mean+bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_std; flipud(bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_mean-bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_std)];        
        
        %plotting and saving
        if(~isnan(bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).mean))
            figure          
            plot(bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).x_unique,bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).mean,'b');
            hold on
            plot(bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).x_rnd_unique,bdt_figures.(struct_tmp{i})(1).testing(plot_id).ROC(j).rnd_mean,'r');
            patch(tmpx_norm, tmpy_norm, tmpx_norm*0-2*eps, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
            
            if(sum(isnan(tmpy_rnd))<length(tmpy_rnd))
                patch(tmpx_rnd, tmpy_rnd, tmpx_rnd*0-2*eps, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2)      
                auc_rnd_text = sprintf('AUC of Random=%g +/- %g',round(100*mean(auc_rnd),2),round(100*std(auc_rnd),2));
                text(0.6,0.05,auc_rnd_text,'EdgeColor','k');
            
            end
            xlabel('False positive rate');
            ylabel('True positive rate');
            auc_text = sprintf('AUC of Training=%g +/- %g',round(100*mean(auc_norm),2),round(100*std(auc_norm),2));
            text(0.6,0.18,auc_text,'EdgeColor','k');
            ylim([0 1]);
            set(gcf,'Color','w');
            %save_name1 = strrep(save_name,'_',' ');        
            title(sprintf('ROC curve of %s for %s',title_tmp,title_name));
            hold off
            saveas(gcf,save_name,'fig')
            saveas(gcf,save_name,'pdf')
        end
    end
end
