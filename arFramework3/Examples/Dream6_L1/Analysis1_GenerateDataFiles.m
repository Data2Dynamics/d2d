%% Selection of experimental conditions
nSimu = 512; % Number of test cases (500 successful - 12 failed)

ratio_condi = .5; % On average, half of conditions observed.

targets = 1:6;
perturbs = 1:3; % 1 = KO, 2 = siRNA, 3 = OE
ratio_mRNA = 1/3; % On average, 1/3 of measurements are mRNA, remaining are protein

rng(0) % Make results reproducible

header = {'Time','pp1_mrna','pp2_mrna','pp3_mrna','pp4_mrna','pp5_mrna','pp6_mrna', ...
    'p1','p2','p3','p4','p5','p6', ...
    'gene1_ko','gene2_ko','gene3_ko','gene4_ko','gene5_ko','gene6_ko', ...
    'rbs1_ic','rbs2_ic','rbs3_ic','rbs4_ic','rbs5_ic','rbs6_ic', ...
    'sirna1_kd','sirna2_kd','sirna3_kd','sirna4_kd','sirna5_kd','sirna6_kd',};
for i = 1:nSimu
    f = fopen(sprintf('Data/simu_%d.csv',i),'w');
    fprintf(f,'%s,',header{1:end-1});
    fprintf(f,'%s\n',header{end});
    fclose(f);
    
    tab = zeros(2,length(header));
    tab(2,1) = 20;
    for j = 1:length(targets)
        for k = 1:length(perturbs)
            if rand < ratio_condi
                % Implement perturbation
                tab = [tab; tab(1:2,:)];
                tab(end-1:end,13+(j-1)*length(perturbs)+k) = 1;
            end
        end
    end
    dlmwrite(sprintf('Data/simu_%d.csv',i),tab,'-append');
    copyfile(sprintf('Data/simu_%d.csv',i),sprintf('Data/simu_%d_l1.csv',i))
    
    copyfile('Data/all_model1.def',sprintf('Data/simu_%d.def',i))
    copyfile('Data/all_model1_l1.def',sprintf('Data/simu_%d_l1.def',i))
end

%% Compile
nSimu = 512;
close all
for i = 1:nSimu
    arInit
    arLoadModel('model1_dream6');
    arLoadData(sprintf('simu_%d',i), 1, 'csv',false);
    arLoadData(sprintf('simu_%d_l1',i), 1, 'csv',false);
    arCompileAll
    ar.config.fiterrors = 0;
    loadGoldStandardParameter(ar.model(1).name) % no priors
    ar.qFit(:) = 1;
    relto = arPrint('relto');
    ar.p(relto) = 0;
    ar.lb(relto(~cellfun(@isempty,strfind(ar.pLabel(relto),'_h')))) = -log10(8);
    ar.ub(relto(~cellfun(@isempty,strfind(ar.pLabel(relto),'_h')))) = log10(8);
    
    arSave(sprintf('simu_%d',i))
end