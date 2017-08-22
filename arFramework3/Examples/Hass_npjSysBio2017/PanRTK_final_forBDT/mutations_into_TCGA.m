% Small script to insert PI3K/RAS mutation status into TCGA matrices

indications = unique(indications_select);

for i=1:length(indications)
    tmp = find(strcmp(indications_select,indications{i})); 
    for j=1:length(tmp)
       bdt.(['TCGA_' indications{i}])((j-1)*bdt.nr_reps+1:j*bdt.nr_reps,bdt.nr_model-1:bdt.nr_model) = repmat(data_select(21:22,tmp(j))',bdt.nr_reps,1);
    end
end