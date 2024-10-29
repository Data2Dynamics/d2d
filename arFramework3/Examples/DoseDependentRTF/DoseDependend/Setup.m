%% This Skript will use the Dose Dependent RTF model to fit all doses at once for the wildtype condition

addpath ../../ % add data2dynamics to path
addpath ../Functions % add functions to path

datafile = 'ELISA-Nigericin-formatted-WT';
model = 'RTF_DoseDependent_WT';
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
arInit; % Initialize d2d
arLoadModel(model); % Load model
arLoadData(datafile); % Load data
arCompileAll; % Compile model
ar.info.datafile = datafile;

%% Set the bounds to the values proposed in the paper
y = [];
for m=1:length(ar.model.data)
    for i=1:length(ar.model.data(m).condition)
        %if ar.model.data(m).condition(i).parameter == "isKO"
          %  if ~str2num(ar.model.data(m).condition(i).value) % only WT
                y = [y, ar.model.data(m).yExp'];
         %   end
        %end
    end
end
t = ar.model.data.tExp;
t = t*10/T;
Delta_t = NaN(length(t)-1,1);
for i=1:(length(t)-1)
    Delta_t(i) = t(i+1)-t(i);
end
% y only for zero dosis
y_0 = [];
for m=1:length(ar.model.data)
    for i=1:length(ar.model.data(m).condition)
        %[y_dose,~] = ar.model.data(m).condition.value
        y_dose = ar.model.data(m).condition.value;
        if y_dose == "1e-05" || y_dose == "0"
            % if ar.model.data(m).condition(i).parameter == "isKO"
            %     if ~str2num(ar.model.data(m).condition(i).value) % only WT
                    y_0 = [y_0, ar.model.data(m).yExp'];
            %     end
            % end
        end
    end
end
doses = doses_num;

% Maximum value of the amplitude of the sustained response
ind = arPrint('M_A');       
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

% Maximum value of the decay rate of the sustained response
ind = arPrint('M_alpha');       
lb  = 1/(2*(max(t)-min(t)));
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

% Maximum value of the time shift
ind = arPrint('M_tau');       
lb = -(max(t)-min(t))/5;
ub = (max(t)-min(t))/2;
p  = -(max(t)-min(t))/10;
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

% Hill coefficient
ind = arPrint('h_');       
%lb = 0.1;
lb = 1;
ub = 10;
p  = 0.5*lb+0.5*ub;
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

% Half maximum quantity
ind = arPrint('K_');       
lb = min(doses(doses>0))/10;
ub = max(doses)*10;
p  = 0.5*lb+0.5*ub;
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

ind = arPrint('sd');       
ub = max(10^-10, std(y));
lb = min(ub, max(10^-10, (max(y) - min(y))/(10^4)));
p  = 0.5*lb+0.5*ub;
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
arSave(datafile,0,0) % save the results

%% Profile Likelihood Estimation (PLE) 
arPLEInit % initialize ple
ple % perform ple
pleSmooth % smooth likelihood profiles
arSave(datafile,0,0) % save the results
