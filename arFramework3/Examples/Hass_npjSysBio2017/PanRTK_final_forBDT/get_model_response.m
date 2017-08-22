%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   crawl_viability                         %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to obtain and normalize fold changes / quasi steady-state / AUC of model components

function get_model_response(whichStruct, whichCell, norm)
global ar
global Matrices
if(~exist('whichStruct','var') || isempty(whichStruct) || strcmp(whichStruct,'ALL') || (isnumeric(whichStruct) && whichStruct==0))
    whichStruct=fetch_outputs;
elseif(ischar(whichStruct))
    which_Struct = whichStruct;
    whichStruct = cell(1);
    whichStruct{1} = which_Struct;
elseif(~ischar(whichStruct) && ~iscell(whichStruct))    
    error('please specifiy the struct to obtain model response, or 0 for ALL');
end

if(~exist('norm','var') || isempty(norm))
   norm = 1;
end

if(~exist('whichCell','var') || isempty(whichCell))
   norm = 1;
end

xlabels = {'control','HRG','HGF','IGF1','HRG/IGF1','EGF','HRG/EGF'};
ylabels = {'media'};

for h=1:length(whichStruct)
    if(~isfield(ar,'AB'))
        ar.AB=[];
    end
    if(~isfield(Matrices,'AB'))
        Matrices.AB=[];
    end
    if(~isfield(ar.AB,whichStruct{h}))
        ar.AB.(whichStruct{h}) = NaN(length(ylabels),length(xlabels),1);
        Matrices.AB.(whichStruct{h}) = NaN(length(ylabels),length(xlabels),1);
    end
    if(~isempty(strfind(whichStruct{h},'qSS')) || ~isempty(strfind(whichStruct{h},'HetHomo')) || ~isempty(strfind(whichStruct{h},'FoldC')))
        norm=0;
        which_course = whichStruct{h};
    else
        which_course = [whichStruct{h} '_course'];
    end
    
    which_ctrlNorm = [whichStruct{h} '_Ctrl_norm'];
    
    ar.AB.(whichStruct{h})(:,:,whichCell)=NaN;
    Matrices.AB.(whichStruct{h})(:,:,whichCell)=NaN;
    for i=1:length(ar.ABplots)
        if(~isfield(ar.ABplots(i),which_course) || sum(sum(isnan(ar.ABplots(i).(which_course))))==length(ar.ABplots(i).(which_course)))
            continue;
        end
        
        celllines_indices = strmatch(ar.ABplots(i).celllines(1),ar.ABplots(i).celllines_names);
        for j=1:length(celllines_indices)
            
            col_tmp = [];
            ks = sum(~cellfun(@isempty,(strfind(ar.model(1).pcond,'_level'))));
            for k=1:ks
                if(str2double(ar.model(1).data(celllines_indices(j)).condition(k).value))
                    if(isempty(col_tmp))
                        col_tmp = strrep(ar.model(1).data(celllines_indices(j)).condition(k).parameter,'_level','');
                    elseif(~strcmp(col_tmp,'HRG'))
                        col_tmp = [strrep(ar.model(1).data(celllines_indices(j)).condition(k).parameter,'_level','') '/' col_tmp];
                    else
                        col_tmp = [col_tmp '/' strrep(ar.model(1).data(celllines_indices(j)).condition(k).parameter,'_level','')];
                    end
                end
            end
            if(isempty(col_tmp))
                col_id = 1;
            else
                col_id = find(ismember(xlabels,col_tmp)==1);
            end
            if(i==1)
                if(norm)
                    ar.AB.(whichStruct{h})(1,col_id,whichCell)=ar.ABplots(i).(which_course)(1,celllines_indices(j))/ar.ABplots(i).(which_ctrlNorm)(1,1);
                    Matrices.AB.(whichStruct{h})(1,col_id,whichCell)=ar.AB.(whichStruct{h})(1,col_id,whichCell);
                else
                    ar.AB.(whichStruct{h})(1,col_id,whichCell)=ar.ABplots(i).(which_course)(1,celllines_indices(j));
                    Matrices.AB.(whichStruct{h})(1,col_id,whichCell)=ar.AB.(whichStruct{h})(1,col_id,whichCell);
                end
            end            
            
        end
    end
end
end