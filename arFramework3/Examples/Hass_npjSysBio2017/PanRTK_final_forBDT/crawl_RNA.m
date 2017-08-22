% Small script to crawl CCLE database for certain receptors/ligand values
% of a list of cell lines

CellLines = RSEM.Properties.VariableNames;
Protein_Names = RSEM{:,1};
Rows = {'EGFR';'ERBB2';'ERBB3';'IGF1R';'ERBB4';'MET';'AKT1';'MAPK1';'HGF';'EGF';'HGF';'HRG';'IGF1';'TGFA';'BTC';'HBEGF'};
tmp = struct();

for i = 1: length(Rows)
    tmp.([Rows{i} '_nr']) = find(strcmp(Rows{i},Protein_Names));
end

%Get complementary cells
ucells_not = find(~ismember(CellLines,ucells2));
ucells_not(1) = [];
nr_not = NaN(length(Rows),length(ucells_not));
for j = 1:length(Rows)
    nr_not(j,:) = RSEM{tmp.([Rows{j} '_nr']),ucells_not};
end

%%
cell_tags = cell(1,length(ucells));
nr_cells = NaN(length(Rows),length(ucells));
ucells2 = ucells;

for i=1:length(ucells)
   tmpcell = strrep(ucells{i},'-','_');
   tmpcell = strrep(tmpcell,'MDA_','MDA_MB_');
   
   tmpcell = strrep(tmpcell,' ','_');
   cellID=regexp(tmpcell,'[a-zA-Z]_?[0-9]{1}');
   tmpcell = strrep(tmpcell,'_','_?');
   if(isempty(cellID))
       tmp1 = sprintf('^.*%s.*$',tmpcell);
   else
    tmp1 = sprintf('^.*%s_?%s.*$',tmpcell(1:cellID),tmpcell(cellID+1:end));
   end
   col_tmp =  find(~cellfun(@isempty,regexpi(CellLines,tmp1)));
   if(isempty(col_tmp))
      sprintf('could not find %s \n',tmp1)
      continue;
   end
   
   if(length(col_tmp)>1)
       sprintf('ohoh');
       for j=1:length(col_tmp)
          CellLines{col_tmp(j)} 
       end
   end
   ucells2{i} = CellLines{col_tmp} ;
    if(strcmp(ucells{i},'H23') || strcmp(ucells{i},'SW48'))
        col_tmp = col_tmp(2);
    end
    if(strcmp(ucells{i},'SNU5'))
        col_tmp = col_tmp(1);
    end
    cell_tags{i} = ucells{i};
    if(length(col_tmp)==1)
        for j = 1:length(Rows)
            nr_cells(j,i) = RSEM{tmp.([Rows{j} '_nr']),col_tmp};
        end
       nr_cells(:,i) = [RSEM{EGFR_nr,col_tmp};RSEM{ERBB2_nr,col_tmp};...
           RSEM{ERBB3_nr,col_tmp};RSEM{IGF1R_nr,col_tmp};RSEM{ERBB4_nr,col_tmp};RSEM{MET_nr,col_tmp};...
           RSEM{HGF_nr,col_tmp};RSEM{EGF_nr,col_tmp};RSEM{HRG_nr,col_tmp};RSEM{IGF_nr,col_tmp}...
           ;RSEM{TGFA_nr,col_tmp};RSEM{BTC_nr,col_tmp};RSEM{HBEGF_nr,col_tmp}];
    end
   
end
empty_cells = find(isnan(nr_cells(1,:)));
ucells2(empty_cells) = [];
nr_cells(:,empty_cells) = [];

