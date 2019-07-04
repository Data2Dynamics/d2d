%% Simulate data at suggested TimePoints
% Load Realistic model
% TimePoints have to be simulated before RealisticSimu with RealisticDesign_D2D
% Or load a realistic model where you just want to change plog or p2

function RealisticSimu(plog,p2)

global ar

if ~exist('plog','var')
    plog = false;
end
if ~exist('p2','var')
    p2 = false;
end
name = 'Realistic';
if p2
    name = [name '_p2'];
end

Magn = ar.Magn;
%init = ar.init;
%yinit = ar.yinit;
tT = ar.model.data.tExp;
if size(tT,2)<1
    tT = xlsread('RealisticDesign/TimePoints.xls');
    if size(tT,2)<1
        error('Something went wrong during RealisticDesign. Sure you loaded a Realistic model before ?!')
    end
end 

%% if Realistic Simu is called for a second time, but this time with log parameters, the same RealisticDesign is used for logarithmic Realisic simulation

if isempty(dir(['Results/' name '_Log'])) && plog
    WriteLogModelDef(name,p2);
    tT = xlsread('RealisticDesign/TimePoints.xls');
elseif ~isempty(dir(['Results/' name '_Log'])) && plog
    warning(['RealisticSimu.m: Logarithmic model already exists. Doesn''t have to be logged.\n'])
    return
end
if plog
    name = [name '_Log'];
end

%% Compile
% has to be! Because error model of sd * Obs and parameters have to be implemented before
% arSimuData
fprintf('Compile new Observables and Parameters.')
arInit
ar.config.checkForNegFluxes = false;

if exist(['Models/' name '.def'],'file')
    arLoadModel(name)
else
    error(['RealisticSimu.m: Couldn''t find ' name '.def']);
end
arLoadData('RealisticData');
arCompileAll(true);
%arSave('Compiled')

% Set tolerances and offset
%ar.config.atol = 10.^(min(Magn))/1000;
%ar.config.rtol = 10.^(min(Magn))/1000;
ar.Magn = Magn;
% ar.init = init;
% ar.yinit = yinit;
ar.config.maxsteps = 10000;
ar.qFit(:)=1;
for i=1:length(ar.p)
    if strncmp(ar.pLabel{i},'offset_',6)
            for j=1:length(ar.model.data.y)
                if strcmp(ar.pLabel{i}(8:end),ar.model.data.y{j}(1:end-4))
                    ar.p(i) = Magn(j)-2;
                    ar.lb(i) = ar.p(i)-2;
                    ar.ub(i) = ar.p(i)+2;
                end
            end
    end
    if strncmp(ar.pLabel{i},'init_',5)
        if plog
            ar.lb(i) = ar.p(i)-1;
            ar.ub(i) = ar.p(i)+1;
        else
            ar.lb(i) = ar.p(i)/10;
            ar.ub(i) = ar.p(i)*10;
        end
    end
%     if strncmp(ar.pLabel{i},'init_',5)
%        for j=1:length(ar.model.x)
%             if strcmp(ar.pLabel{i}(6:end),ar.model.x{j})
%                 if plog
%                     ar.p(i) = log10(init(j));
%                     ar.lb(i) = ar.p(i)-2;
%                     ar.ub(i) = ar.p(i)+2;
%                     ar.qLog10(i)=1;
%                 else
%                     ar.p(i) = init(j);
%                     ar.lb(i) = ar.p(i)/100;
%                     ar.ub(i) = ar.p(i)*100;
%                     ar.qLog10(i)=0;
%                 end
%             end
%        end
%     end
    if strncmp(ar.pLabel{i},'scale_',6)     % noch bei -3,3 bei bisher compiled modellen
        ar.lb(i) = -2;
        ar.p(i) = 0;
        ar.ub(i) = 2;
    end
    if strncmp(ar.pLabel{i},'sd_',3)        
        %arLink
        %arSimu(true,true)
        if strncmp(ar.pLabel{i},'sd_rel',6)
            ar.p(i) = -1;
            ar.lb(i) = -2;
            ar.ub(i) = 0;
        else
            for j=1:length(ar.model.data.y)
                if strcmp(ar.pLabel{i}(4:end),ar.model.data.y{j})
                    ar.p(i) = Magn(j)-1;
                    ar.lb(i) = ar.p(i)-1;
                    ar.ub(i) = ar.p(i)+1;
                end
            end
        end
    end
end


arPlot
ar.config.fiterrors = 1; 

%arSave('Pars')
% tT = xlsread('RealisticDesign/TimePoints.xls');
y = nan(size(tT,1),size(tT,2));
for i = 1:size(tT,2) 
    ar.model.data.tExp = tT(~isnan(tT(:,i)),i);
    ar.model.data.yExp = [];
    arSimuData;
    arPlot
%       print('AppendInSimu','-dpdf')
%       append_pdf('RealisticDesign/Simu.pdf','Simu.pdf')
%       print -dpsc -append RealisticDesign/Simu.ps
%   Just remember Observable i with the specific time points tT(i)
    y(1:size(ar.model.data.tExp,1),i) = ar.model.data.yExp(:,i);
end
% system('ps2pdf RealisticDesign/Simu.ps');

WriteDataTableRealistic(tT, y, ar.model.data.y);

%% Set Init Parameters to first value of simulated data
% for i=1:length(ar.p)
%     if strncmp(ar.pLabel{i},'init_',5)
%        for j=1:length(ar.model.data.y)
%             if strcmp(ar.pLabel{i}(6:end),ar.model.data.y{j})
%                 if plog
%                     ar.p(i) = log10(ar.model.data.yExp(1,j));
%                     ar.lb(i) = ar.p(i)-2;
%                     ar.ub(i) = ar.p(i)+2;
%                     ar.qLog10(i)=1;
%                 else
%                     ar.p(i) = ar.model.data.yExp(1,j);
%                     ar.lb(i) = ar.p(i)/100;
%                     ar.ub(i) = ar.p(i)*100;
%                     ar.qLog10(i)=0;
%                 end
%             end
%        end
%     end
% end

ar.config.ploterrors=2;
arLink;
arSimu(true,true);
arPlot;
if exist(['Results\' name],'dir')
    rmdir(['Results\' name],'s');
end
arSave(name)
[~,ws]=fileparts(ar.config.savepath);
movefile(['Results/' ws],['Results\' name]);
fprintf(['Realistic simulation saved to ./Results/' name '/workspace.mat \n']);

