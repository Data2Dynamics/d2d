clear all;
close all;
clc;

%% Load Model
arInit
arLoadModel('AktPathwayFujita'); model_name = 'AktPathwayFujita';

%% Load Data
dsl = {}; dls_filename = {};
f = dir(['./Data/*.def']); file_list = {f.name}'; file_list = natsortfiles(file_list);

for n = 1:6%length(file_list)
    dls_filename{end+1} = file_list{n}(1:end-4);
    dsl{end+1}.name = {file_list{n}(1:end-4)};
end

% Loop: Datasets and Conditions
warning off
for i = 1:length(dsl)
    for j = 1:length(dsl{i}.name)
        arLoadData(dsl{i}.name{j},1,'csv');
    end
end
warning on

% Compile model
warning off
arCompileAll;
warning on

% Show condition
arShowDataConditionStructure

%% Constraint parameters
ar.lb = ar.lb.*0-8;
ar.ub = ar.ub.*0+8;

%% Numerical settings
ar.config.atol = 1e-6;
ar.config.rtol = 1e-6;
ar.config.maxsteps = 1e4;

ar.config.optim = optimset(...
    ar.config.optim,'TolFun',10^-8,'PrecondBandWidth',inf,...
    'display','iter');
ar.config.optim.MaxIter = 2000;

%% Set pre-equilibration
% ic_MKN1_FM   = arFindCondition(ar,'experimentaldata'  ,'EGF_conc_step'  ,0,'EGF_conc_impulse'  ,0,'EGF_conc_ramp',0,'pulse_time',0,'ramp_time',0);
% 
% if ~isempty(ic_MKN1_FM)
%     arSteadyState(1,ic_MKN1_FM  ,  arFindCondition(ar,'experimentaldata'  ,'EGF_conc_step'  ,0,'EGF_conc_impulse'  ,0,'EGF_conc_ramp',0,'pulse_time',0,'ramp_time',0));
% end

arSimu(true,true,true);

%% Fitting
n_fit = 50;
if n_fit >= 2
    arFitLHS(n_fit);
end
save('results.mat')
theta_optim = ar.p;

%% Check paper parameters
arFits(theta_optim)

theta_true = [1.06386e-4 4.3309e-2 6.81902e4 3.54317 6.00588e1 1.81735e-4 4.98862e4 6.73816e-3 4.0749e-2 1.92391e-2 9.97194e-2 1.5543e-5 5.17473e-3 3.05684e-2 3.27962e-2 2.10189e-6 5.1794e-15 1.21498e-3 1.13102e-3 0.1 0.1 0.1];
theta_true = log10(theta_true);
arFits(theta_true)

%% Visualization options
fbound = 1.5;
n_sigma = 2;
d_min = 1e-4;
fs = 8;
w_std = 0.2;
% ar.p = theta_optim;
arSimu(true,true,true)

y_name ={'pEGFR','pAkt','pS6'};
close all
% Loop: Datasets
% Initialization
leg_h = nan(length(dsl),1);
leg_n = {};
for d = 1:length(dsl)
%     figure('name',[strrep(dls_filename{d},'_',' ') ' - fit']);

    % Compute dataset index
    i = 0;
    for delta = 1:d-1
        i = i + length(dsl{delta}.name);
    end
    
    % Colormap
    col_sim = colormap(jet(256));
    col_data = colormap(jet(256));
    col_sim = col_sim(1:42:end,:);
    col_data = col_data(1:42:end,:);
    
    % Datatype
    switch ar.model.plot(i+1).doseresponse
        case 0 % Kinetic
            % Loop: Conditions of dataset d
            for c = 1:length(dsl{d}.name)
                % Data and condition index
                Id = ar.model.plot(i+c).dLink;
                Ic = ar.model.data(Id).cLink;

                for y_ind = 1:length(y_name)
                    subplot(2,3,y_ind);

                    % Collection of data
                    y = [];
                    t_sim = ar.model.condition(Ic).tFine;
                    for y_ind_Id = 1:length(ar.model.z)
                        if strfind(ar.model.z{y_ind_Id},y_name{y_ind}) == 1
                            y_sim = ar.model.condition(Ic).zFineSimu(:,y_ind_Id);
                        end
                    end
                    y_exp = [];y_exp_std = [];
                    cy_ind_Id = 1;
                    for y_ind_Id = 1:length(ar.model.data(Id).yNames)
                        if strfind(ar.model.data(Id).yNames{y_ind_Id},y_name{y_ind}) == 1
                            sc_y = 1;
%                             sc_y = arGetPars(['s_' ar.model.data(Id).yNames{y_ind_Id}],0);
                            y(:,cy_ind_Id) = ar.model.data(Id).yExpSimu(:,y_ind_Id)/sc_y;
                            y_exp(:,cy_ind_Id) = ar.model.data(Id).yExp(:,y_ind_Id)'/sc_y;
                            y_exp_std(:,cy_ind_Id) = ar.model.data(Id).yExpStd(:,y_ind_Id);
                            cy_ind_Id = cy_ind_Id + 1;
                        end
                    end
                        
                    % Experimental data
                    if c == 1
                        w_std_t = diff(ar.model.data(Id).tLim)/50;
                    end
                    
                    if ~isempty(y_exp)
                        leg_h(d) = plot(ar.model.data(Id).tExp,nanmean(y_exp,2),'o--','color',col_data(d,:),'linewidth',1); hold on;
                        for k = 1:length(ar.model.data(Id).tExp)
                            errorbar(ar.model.data(Id).tExp(k),y_exp(k,:),y_exp_std(k,:),'-','color',col_data(d,:),'linewidth',1); hold on;
%                             plot(ar.model.data(Id).tExp(k)+w_std_t*[-1,1],nanmean(y_exp(k,:))-n_sigma*nanstd(y_exp(k,:))*[1,1],'-','color',col_data(c,:),'linewidth',2); hold on;
%                             plot(ar.model.data(Id).tExp(k)+w_std_t*[-1,1],nanmean(y_exp(k,:))+n_sigma*nanstd(y_exp(k,:))*[1,1],'-','color',col_data(c,:),'linewidth',2); hold on;
                        end
%                         leg_n{d} = [ar.model.data(d).condition(1).value ' - Data'];
                    end
                    
                    % Simulation
                    leg_h(d) = plot(t_sim,y_sim,'-','color',col_sim(d,:),'linewidth',1); hold on;
                    leg_n{d} = ['EGF ' ar.model.data(d).condition(1).value 'ng/ml'];
                    
                    % Label
                    if c == 1
                        xlabel('time [s]')
                        ylabel(y_name{y_ind})
                        set(gca,'fontsize',fs);
                    end
                    
                    % Limits
                    if c == length(dsl{d}.name)
                        cxlim = ar.model.data(Id).tLim;
                        cylim = get(gca,'ylim'); cylim(1) = 0;
                        set(gca,'ylim',[0 1.25],'xlim',[0 3600]);
                    end
                end
            end
    end        
end
% Legend & title
legend(leg_h(find(leg_h~=0)),leg_n{find(leg_h~=0)});

    % Save figure
%     set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 12])
%     print('-depsc2','-r1200',['./Figures/' model_name '__' dls_filename{d} '__fit']);
