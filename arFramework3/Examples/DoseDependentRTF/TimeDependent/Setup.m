%% This Skript will use the Single Dose RTF model to fit the data seperately for each dose for the wildtype condition

addpath ../../ % add data2dynamics to path
addpath ../Functions % add functions to path

datafile = 'ELISA-Nigericin-formatted' ;
model = 'RTF_TimeDependent';
T = arSetTRange(model, datafile, 'T', 0, 'data'); % Find range of the measurement times T and write it into model def file

%% Get some information about the doses from the data file
[d,s] = xlsread(['Data/',datafile,'.xlsx']);
s = s(1,:);
index_dose = find(contains(s,'dose'));
doses_len = length(unique(d(:,index_dose)));
doses_num = unique(d(:,index_dose));
doses_str = num2str(doses_num);
filesDoses = {};
for i = 1:doses_len;
    i_str = strtrim(doses_str(i,:));
    filesDoses{i} = ['Doses/',datafile,'-Doses/',datafile,'-Dosis-', i_str];
end

%% Beginn d2d modelling
for k=1:doses_len
    arInit; % Initialize d2d
    arLoadModel(model); % Load model
    arLoadData(filesDoses{k}); % Load data
    arCompileAll; % Compile model

    ar.doses = unique(d(:,index_dose)); 
    ar = orderfields(ar);

    %% Set the bounds to the values proposed in the paper

    y = ar.model.data.yExp(~isnan(ar.model.data.yExp)); % get the experimental data
    t = unique(ar.model.data.tExp*10/T); % get the experimental time points
    Delta_t = NaN(length(t)-1,1); 
    for i=1:(length(t)-1)
        Delta_t(i) = t(i+1)-t(i); % calculate the time differences
    end

    % Amplitude of the sustained response
    ind = arPrint('A');       
    % lb  = -2*(max(y)-min(y)); 
    lb  = 0; 
    ub  = 2*(max(y)-min(y)); 
    p   = 0.1*lb+0.9*ub;
    if lb<=0
        ar.qLog10(ind) = 0;
        ar.ub(ind) = ub;
        ar.lb(ind) = lb;
        ar.p(ind) = p;   
    else
        ar.ub(ind) = log10(ub);
        ar.lb(ind) = log10(lb);
        ar.p(ind) = log10(p);   
    end

    % Decay rate of the sustained response
     ind = arPrint('alpha');       
    lb  = 1/(max(Delta_t)*2);
    ub  = 2/min(Delta_t);
    p   = 0.5*lb+0.5*ub;
    if lb<=0
        ar.qLog10(ind) = 0;
        ar.ub(ind) = ub;
        ar.lb(ind) = lb;
        ar.p(ind) = p;   
    else
        ar.ub(ind) = log10(ub);
        ar.lb(ind) = log10(lb);
        ar.p(ind) = log10(p);   
    end

    % Offset b
    ind = arPrint('b');       
    lb  = min(y);
    ub  = max(y);
    p   = 0.5*lb+0.5*ub;
    if lb<=0
        ar.qLog10(ind) = 0;
        ar.ub(ind) = ub;
        ar.lb(ind) = lb;
        ar.p(ind) = p;   
    else
        ar.ub(ind) = log10(ub);
        ar.lb(ind) = log10(lb);
        ar.p(ind) = log10(p);   
    end

    % Time shift
    ind = arPrint('tau');       
    lb  = -(max(t)-min(t))/5;
    ub  = (max(t)-min(t))/2;
    p   = -(max(t)-min(t))/10;
    if lb<=0
        ar.qLog10(ind) = 0;
        ar.ub(ind) = ub;
        ar.lb(ind) = lb;
        ar.p(ind) = p;   
    else
        ar.ub(ind) = log10(ub);
        ar.lb(ind) = log10(lb);
        ar.p(ind) = log10(p);   
    end

    arPrint % print the parameters with the new bounds
    
    %% Fit the model to the data
    arFitLHS(100) % multistart optimization
    arWaterfallPlot % waterfall plot
    arFit % fit model to data
    arPlot % plot model and data
    arPrint % print estimated parameters  
    arSave([datafile,'_dose_',num2str(k)],0,0) % save the results

    %% Profile Likelihood Estimation (PLE)
    arPLEInit % initialize ple
    ple % perform ple
    pleSmooth % smooth likelihood profiles
    arSave([datafile,'_dose_',num2str(k)],0,0) % save the results
end
