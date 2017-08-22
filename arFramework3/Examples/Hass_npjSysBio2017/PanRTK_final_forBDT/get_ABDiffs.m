%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   get_ABDiffs                             %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get AUC / quasi SS/ fold-change for receptors and downstream components
%  for all cell lines and ligand conditions

%Ligands: 
%HRG       5nN
%HGF       1nM
%EGF        5nM
%FGF2      5nM
%IGF1       50nM 
%IGF2       50nM
%Insulin   50nM

%add tExp, t at ~4000 min, tLim(2) = 4000

function get_ABDiffs(tEnd)
global ar

if(~exist('tEnd','var') || isempty(tEnd))
    tEnd=1440;
end

ar.ABplots.input=0;

%Extend time frame to provided read-out time
reset_tLim = ar.model(1).data(1).tLim;
if(ar.model(1).data(1).tExp(end)~=tEnd)
   
    for i=1:length(ar.model(1).data)
        ar.model(1).data(i).tLim(2) = tEnd;
        ar.model(1).tLim(2) = tEnd;
        ar.model(1).data(i).tExp=[ar.model(1).data(i).tExp;tEnd];
        ar.model(1).data(i).yExp=[ar.model(1).data(i).yExp; NaN(1,size(ar.model(1).data(i).yExp,2))];
        ar.model(1).data(i).yExpStd=[ar.model(1).data(i).yExpStd; NaN(1,size(ar.model(1).data(i).yExp,2))];
    end
    arLink
end

%Collect indices of targets
pAKT_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pAKT')));
AKT_index = find(~cellfun(@isempty,regexp(ar.model(1).xNames,'^(AKT)*')));
pERK_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pERK')));
ERK_index = find(~cellfun(@isempty,regexp(ar.model(1).xNames,'^(ERK)*')));
pS6_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pS6')));
pS6_index = pS6_index(2);

EGFR_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'EGFR'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'EGFR'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB1'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'MetEGFR')))];

pEGFR_hom_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pEGFR')))];
pEGFR_12_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB12')));
pEGFR_13_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB13')));
pEGFR_Met_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetEGFR')));
pEGFR_het_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB1'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetEGFR')))];

ErbB3_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB3'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB3')) & cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB32'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB13'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'MetErbB3')))];

pErbB3_hom_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB3')) & cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB32')));
pErbB3_32_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB32')));
pErbB3_13_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB13')));
pErbB3_Met_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetErbB3')));
pErbB3_het_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB32'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB13'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetErbB3')))];

ErbB2_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB2'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB2'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB12'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'ErbB32')))];

pErbB2_hom_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB2')));
pErbB2_12_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB12')));
pErbB2_32_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB32')));
pErbB2_het_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB12'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pErbB32')))];

Met_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'Met')))];

pMet_hom_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMet')) & cellfun(@isempty,strfind(ar.model(1).xNames,'pMetEGFR')) & cellfun(@isempty,strfind(ar.model(1).xNames,'pMetErbB3')))];
pMet_EGFR_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetEGFR')));
pMet_ErbB3_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetErbB3')));
pMet_het_index = [find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetEGFR'))) find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pMetErbB3')))];

pIGF1R_hom_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'pIGF1R')));
IGF1R_index = find(~cellfun(@isempty,strfind(ar.model(1).xNames,'IGF1R')));

pEGFR_index = [pEGFR_hom_index pEGFR_het_index];
pErbB3_index = [pErbB3_hom_index pErbB3_het_index];
pErbB2_index = [pErbB2_hom_index pErbB2_het_index];
pMet_index = [pMet_hom_index pMet_het_index];

list_outputs = fetch_outputs;
input = 0;

ar.ABplots.celllines_names=cell(length(ar.model(1).data),1);
for i=1:length(ar.model(1).data)
    ar.ABplots.celllines_names{i} = ar.model(1).data(i).name;
end
ar.ABplots.celllines = unique(ar.ABplots.celllines_names);
celllines_indices = strmatch(ar.ABplots.celllines(1),ar.ABplots.celllines_names);
    
for i=1:length(list_outputs)
    if(~isempty(strfind(list_outputs{i},'qSS')) || ~isempty(strfind(list_outputs{i},'HetHomo')) || ~isempty(strfind(list_outputs{i},'FoldC')))
        tmp_course = list_outputs{i};
        ar.ABplots.(tmp_course) =  NaN(length(input),length(ar.model(1).condition));
    else
        tmp_course = [list_outputs{i} '_course'];
        tmp_ctrlnorm = [list_outputs{i} '_Ctrl_norm'];
        ar.ABplots.(tmp_course) =  NaN(1,length(ar.model(1).condition));
        ar.ABplots.(tmp_ctrlnorm) = NaN(1,length(ar.ABplots.celllines));        
    end
    
end

for j=1:length(ar.ABplots.input)
    
    norm_index = 1;   
    try
        arSimu(false,true,true)
    catch
        sprintf('Trying another Simu run')
        if(ar.config.atol>1.e-9)
            ar.config.atol=0.5*ar.config.atol;
            ar.config.rtol=0.5*ar.config.rtol;
        end
        try
            arSimu(false,true,true)
        catch
            error('simulation with given expression did not work, try other configurations or check RTK expressions')
        end
    end
    for i=1:length(ar.model(1).condition)
        if(isempty(pMet_hom_index))
          continue; 
        end
        if(strcmp(ar.model(1).data(ar.model(1).condition(i).dLink).condition(1).value,'0') && strcmp(ar.model(1).data(ar.model(1).condition(i).dLink).condition(2).value,'0') && strcmp(ar.model(1).data(ar.model(1).condition(i).dLink).condition(3).value,'0') && strcmp(ar.model(1).data(ar.model(1).condition(i).dLink).condition(4).value,'0') && j==1)                 
           for k=1:length(list_outputs)
               if(~isempty(strfind(list_outputs{k},'qSS')) || ~isempty(strfind(list_outputs{k},'HetHomo')) || ~isempty(strfind(list_outputs{k},'FoldC')))
                    continue
               end
               tmp_ctrlnorm = [list_outputs{k} '_Ctrl_norm'];
               tmp_index = [list_outputs{k} '_index'];          
               ar.ABplots.(tmp_ctrlnorm)(1,norm_index) = nansum(nansum(ar.model(1).condition(i).xFineSimu(:,eval(tmp_index))))/(length(eval(tmp_index)));                

            end
           norm_index=norm_index+1;
        end       

        for k=1:length(list_outputs)           
            if(~isempty(strfind(list_outputs{k},'qSS')) || ~isempty(strfind(list_outputs{k},'HetHomo')) || ~isempty(strfind(list_outputs{k},'FoldC')))
                which_index = strsplit(list_outputs{k},'_');
                tmp_index = [which_index{1} '_index'];
                tmp_course = list_outputs{k};
            else
                tmp_index = [list_outputs{k} '_index'];
                tmp_course = [list_outputs{k} '_course'];
            end
            if(~isempty(strfind(list_outputs{k},'qSS')))
                ar.ABplots.(tmp_course)(j,i) = nansum(ar.model(1).condition(i).xFineSimu(end,eval(tmp_index)),2)/(nansum(ar.model(1).condition(i).xFineSimu(1,eval(tmp_index)),2));
            elseif(~isempty(strfind(list_outputs{k},'FoldC')))
                ar.ABplots.(tmp_course)(j,i) = max(nansum(ar.model(1).condition(i).xFineSimu(:,eval(tmp_index)),2))/(nansum(ar.model(1).condition(i).xFineSimu(1,eval(tmp_index)),2));
            else
                ar.ABplots.(tmp_course)(j,i) = nansum(nansum(ar.model(1).condition(i).xFineSimu(:,eval(tmp_index))))/(length(eval(tmp_index)));
            end
        end        
    end     
end


for i=1:length(ar.model(1).data)
    ar.model(1).data(i).tExp(end)=[];
    ar.model(1).data(i).yExp(end,:) = [];
    ar.model(1).data(i).yExpStd(end,:) = [];
    ar.model(1).data(i).tLim=reset_tLim;
end
ar.model(1).tLim = reset_tLim;

