% Looping over actual conditions

close all
results = dir('Results');
targets = 1:6;
ratio_mRNA = 1/3; % On average, 1/3 of measurements are mRNA, remaining are protein

fcs = [2 5 10];
fcs = [1./fcs fcs];
fcs_h = [2 4];
fcs_h = [1./fcs_h fcs_h];
ratio_fcs = 1/3; % On average, 1/3 of parameters are different between cell types
fcs = log10(fcs);
fcs_h = log10(fcs_h);

linv = 1./10.^[-4 -3 -2 -1 0 1 1.1:.1:2.9 3 4 5 6];

for i = 1:length(results)
    rng(i) % Make results reproducible
    
    if ~isempty(strfind(results(i).name,'_simu_'))
        arLoad(results(i).name)
        
        global arWaitbarGlobal;
        arWaitbarGlobal.showWindow = 0;
        
        relto = arPrint('relto');
        for ip = 1:length(relto)
            if rand < ratio_fcs
                p = ar.p(strcmp(ar.pLabel,ar.pLabel{relto(ip)}(7:end)));
                if ~isempty(strfind(ar.pLabel{relto(ip)},'_h'))
                    % Hill (restrict relto to [1 4])
                    myfcs = fcs_h(p+fcs_h >= 0 & p+fcs_h <= log10(4));
                else
                    % Not Hill
                    myfcs = fcs;
                end
                myper = randperm(length(myfcs));
                ar.p(relto(ip)) = myfcs(myper(1));
            end
        end
        arSimuData(1,1:2,0:.5:20)
        for id = 1:length(ar.model.data)/2
            if rand < ratio_mRNA
                % Observe mRNA (21 data-points)
                ar.model.data(id).yExp(:,7:12) = nan;
                ar.model.data(id+length(ar.model.data)/2).yExp(:,7:12) = nan;
                ar.model.data(id).yExp = ar.model.data(id).yExp(1:2:end,:);
                ar.model.data(id).tExp = ar.model.data(id).tExp(1:2:end);
                ar.model.data(id+length(ar.model.data)/2).yExp = ar.model.data(id+length(ar.model.data)/2).yExp(1:2:end,:);
                ar.model.data(id+length(ar.model.data)/2).tExp = ar.model.data(id+length(ar.model.data)/2).tExp(1:2:end);
            else
                % Observe protein
                mytar = randperm(length(targets));
                % Observe mytar(1:2) (41 data-points)
                ar.model.data(id).yExp(:,[1:6 6+mytar(1:4)]) = nan;
                ar.model.data(id+length(ar.model.data)/2).yExp(:,[1:6 6+mytar(1:4)]) = nan;
            end
        end
        arLink
        
        try
            arFit
            l1pars = relto;
            l1Init(l1pars,l1pars*0,ar.lb(l1pars),ar.ub(l1pars))
            l1Scan(l1pars,linv)
            l1Unpen(l1pars)
            l1SelectOpt

            arSave(sprintf('result_%s',results(i).name(22:end)))
        catch exception
            fprintf('%s\n', exception.message);
            arSave(sprintf('result_%s_failed',results(i).name(22:end)))
        end
        
    end
end