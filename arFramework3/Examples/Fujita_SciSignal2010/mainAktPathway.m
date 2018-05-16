% Main file of the AktPathway example
%
% Demonstrates the use of:
% * getMultiStarts()
% * getParameterProfiles()
%
% This model is taken from the paper "Decoupling of Receptor and Downstream
% Signals in the Akt Pathway by Its Low-Pass Filter Characteristics" , by
% Fujita et. al., in 2010, ScienceSignaling, Vol.3 Issue 132 ra56 (see
% https://www.ncbi.nlm.nih.gov/pubmed/20664065)

% The data used is measurement data visualized in the publication. It was
% extracted from plots using WebPlotDigitizer
%
% This file performs a multistart local optimization based on measured data 
% from the referenced paper, demonstrating the use of getMultiStarts().
% Afterwards, the parameter profiles are computed to validate the choice of
% the model. Finally, the simulated data is plotted against the 
% experimental data to visually validate the parameters that were found.


%% Preliminary
% Clean up

clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

% Seed random number generator
rng(0);

%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see EGFR_model_syms.m


[exdir,~,~] = fileparts(which('mainAktPathway.m'));
try
    amiwrap('AktPathway_model', 'AktPathway_syms', exdir);
catch ME
    warning('There was a problem with the AMICI toolbox (available at https://github.com/ICB-DCM/AMICI), which is needed to run this example file. The original error message was:');
    rethrow(ME);
end

%% Data
% Experimental data is read out and written to an AMICI-data object which 
% is used for the ODE integration

D = getData_Akt_pathway();

%% Generation of the structs and options for PESTO
% The structs and the PestoOptions object, which are necessary for the 
% PESTO routines to work are created and set to convenient values


% Parameters that were found initially by Fujita et al.
theta_true = log10([0.00673816; 0.0407490; 0.0192391;...
    0.0997194; 0.0000155430;0.00517473; 0.0305684; 0.0327962; 0.00000210189;...
    0.00000000000000517940; 0.00121498;0.00113102; 0.000106386;0.000181735; 60.00588; 49886.2]);

% Write the parameters struct
parameters.name = {'EGF+EGFR\_k1' 'EGF+EGFR\_k2' 'EGFR\_phosphorylation\_k1' 'pEGFR\_degradation\_k1' 'pEGFR+Akt\_k1' 'pEGFR+Akt\_k2' 'Akt\_phosphorylation\_k1' 'pAkt\_dephospho\_k1' 'pAkt+S6\_k1' 'pAkt+S6\_k2' 'S6\_phosphorylation\_k1' 'pS6\_dephospho\_k1' 'EGFR\_turnover' 'pEGFR\_scaleFactor' 'pAkt\_scaleFactor' 'pS6\_scaleFactor'} ;
parameters.min=[-5,-25,-6,-5,-15,-20,-15,-7,-18,-55,-8,-7,-13,-25,-2,-2];
parameters.max=[3,1,10,4,1,6,3,5,2,5,2,1,1,0,8,25];
parameters.number = length(parameters.name);


% Set the objective function
objectiveFunction = @(theta) LogLikelihood_AktPathway(theta, D(1:6));


% Set the PESTO-options
optionsMultistart           = PestoOptions();
optionsMultistart.n_starts  = 1000;
optionsMultistart.comp_type = 'sequential';
optionsMultistart.mode      = 'visual';
optionsMultistart.obj_type = 'log-posterior';
optionsMultistart.plot_options.add_points.par = theta_true;
optionsMultistart.plot_options.add_points.logPost = objectiveFunction(theta_true);
optionsMultistart.localOptimizerOptions = optimset('Algorithm','interior-point',...
    'GradObj', 'on',...
    'Display', 'iter', ...
    'MaxIter', 2000,...
    'TolFun', 0,...
    'TolX', 1e-10,...
    'MaxFunEvals', 2000);


%% test Gradient
% activate this section to compute finite difference approximations to the
% gradient to check whether the gradient is computed in the right way

%thet=rand(1,18);
%[g,g_fd_f,g_fd_b,g_fd_c]=testGradient(parameters.MS.par(:,1),objectiveFunction,1e-5);


%% Perform Multistart optimization
% A multi-start local optimization is performed within the bound defined in
% parameters.min and .max in order to infer the unknown parameters from 
% measurement data.

fprintf('\n Perform optimization...');
parameters = getMultiStarts(parameters, objectiveFunction, optionsMultistart);


%% Parameter profiles
% The parameter profiles are computed in order to see whether the
% parameters are identifiable and thus the model is well chosen.

parameter_profiles=getParameterProfiles(parameters,objectiveFunction,optionsMultistart);
%% PLOTTING
% loading the parameters from the optimazation
p = parameters.MS.par(:,1);  

% Plotting of experimental Data and simulated Data in one PLOT  only for
% STEP (log10 transformed)

figure(1)
t = linspace(0,6000,6001);
options = [];
tit = ["pEGFR", "pAkt", "pS6"];
t_new = log10(t+1);

t_points=log10([0,1,2,5,10,15,30,60]*60+1); 
sol_p1= simulate_AktPathway_model(t,p,D(1).condition,[],options);
sol_p2= simulate_AktPathway_model(t,p,D(2).condition,[],options);
sol_p3= simulate_AktPathway_model(t,p,D(3).condition,[],options);
sol_p4= simulate_AktPathway_model(t,p,D(4).condition,[],options);
sol_p5= simulate_AktPathway_model(t,p,D(5).condition,[],options);
sol_p6= simulate_AktPathway_model(t,p,D(6).condition,[],options);
        
    for i = 1:3
        m = ((i-1)*3)+1;
        subplot(3,1,i);

            plot(t_new, sol_p1.y(:,i), 'k',...
            t_new, sol_p2.y(:,i),'b',...
            t_new, sol_p3.y(:,i),'g',...
            t_new, sol_p4.y(:,i),'c',...
            t_new, sol_p5.y(:,i),'y',...
            t_new, sol_p6.y(:,i),'r',...
            t_points, D(1).Y(:,i),'ok',...
            t_points, D(2).Y(:,i),'ob',...
            t_points, D(3).Y(:,i),'og',...
            t_points, D(4).Y(:,i),'oc',...
            t_points, D(5).Y(:,i),'oy',...
            t_points, D(6).Y(:,i),'or','LineWidth',2, 'Markersize', 4);

        switch i
           case 1
                axis([0 log10(6001) 0 1.1]);
                leg = legend('0.1', '0.3','1','3','10','30');
                leg.FontSize =14;
                leg.Units = 'centimeters';
                leg.Orientation = 'horizontal';
                leg.Location = 'North';
                legend boxoff
                title(leg,'EGF (ng/ml) in step stimulation','FontSize',14)
            case 2
                axis([0 log10(6001) 0 1.1]);
            case 3
                axis([0 log10(6001) 0 1.0]);
        end
        set(gca,'FontSize',14) 
        xticks([log10(1) log10(61) log10(601) log10(6001)]);
        xticklabels({'0','1','10','100'});
        ylabel('Phosphorylation(AU)','FontSize',14);
        xlabel('time','FontSize',14);
        title(tit(i),'FontSize',14);
    end

%%
% Plotting of experimental Data and simulated Data in one PLOT for ramp and
% pulse (with log10)
figure(2);
t = linspace(0,6000,6001);
options = [];
tit = ["pEGFR", "pAkt", "pS6"];
t_new = log10(t+1);

for j = 2:3
    switch j 
        case 2
        t_points=log10([0,1,2,5,10,15,30,60]*60+1); 
        sol_p1= simulate_AktPathway_model(t,p,D(1+(j-1)*6).condition,[],options);
        sol_p2= simulate_AktPathway_model(t,p,D(2+(j-1)*6).condition,[],options);
        sol_p3= simulate_AktPathway_model(t,p,D(3+(j-1)*6).condition,[],options);
        sol_p4= simulate_AktPathway_model(t,p,D(4+(j-1)*6).condition,[],options);
        sol_p5= simulate_AktPathway_model(t,p,D(5+(j-1)*6).condition,[],options);
        sol_p6= simulate_AktPathway_model(t,p,D(6+(j-1)*6).condition,[],options);
        case 3
        t_points=log10([0,1,2,5,10,15,20,30,40,50,60]*60+1); 
        sol_p1= simulate_AktPathway_model(t,p,D(1+(j-1)*6).condition,[],options);
        sol_p2= simulate_AktPathway_model(t,p,D(2+(j-1)*6).condition,[],options);
        sol_p3= simulate_AktPathway_model(t,p,D(3+(j-1)*6).condition,[],options);
        sol_p4= simulate_AktPathway_model(t,p,D(4+(j-1)*6).condition,[],options);
    end
        for i = 1:3
            m = (i-1)*2+(j-1);
            subplot(3,2,m);
            switch j
             case  {2}
                plot(t_new, sol_p1.y(:,i), 'k',...
                t_new, sol_p2.y(:,i),'b',...
                t_new, sol_p3.y(:,i),'g',...
                t_new, sol_p4.y(:,i),'c',...
                t_new, sol_p5.y(:,i),'y',...
                t_new, sol_p6.y(:,i),'r',...
                t_points, D(1+(j-1)*6).Y(:,i),'ok',...
                t_points, D(2+(j-1)*6).Y(:,i),'ob',...
                t_points, D(3+(j-1)*6).Y(:,i),'og',...
                t_points, D(4+(j-1)*6).Y(:,i),'oc',...
                t_points, D(5+(j-1)*6).Y(:,i),'oy',...
                t_points, D(6+(j-1)*6).Y(:,i),'or','LineWidth',2, 'Markersize', 4);           
            switch i
                case 3
                                axis([0 log10(6001) 0 1.0]);
                case {1,2}
                                axis([0 log10(6001) 0 1.1]);
            end
             case  3 
                plot(t_new, sol_p1.y(:,i), 'b',...
                t_new, sol_p2.y(:,i),'g',...
                t_new, sol_p3.y(:,i),'y',...
                t_new, sol_p4.y(:,i),'r',...
                t_points, D(1+(j-1)*6).Y(:,i),'ob',...
                t_points, D(2+(j-1)*6).Y(:,i),'og',...
                t_points, D(3+(j-1)*6).Y(:,i),'oy',...
                t_points, D(4+(j-1)*6).Y(:,i),'or','LineWidth',2, 'Markersize', 4);
                axis([0 log10(6001) 0 1.0]);
            end
            switch j
                case 2
                    switch i
                        case 1
                        leg = legend('0.1', '0.3','1','3','10','30');
                        leg.FontSize =14;
                        leg.Units = 'centimeters';
                        leg.Orientation = 'horizontal';
                        leg.Location = 'North';
                        title(leg,'EGF (ng/ml) in pulse stimulation','FontSize',14)
                        legend boxoff 
                        end
                case 3
                    switch i
                        case 1
                            leg = legend('0.03', '0.3','3','30');
                            leg.FontSize =14;
                            leg.Units = 'centimeters';
                            leg.Orientation = 'horizontal';
                            leg.Location = 'North'; 
                            title(leg,'EGF (ng/ml) in ramp stimulation','FontSize',14)
                            legend boxoff 
                        end
            end
            set(gca,'FontSize',14) 
            xticks([log10(1) log10(61) log10(601) log10(6001)]);
            xticklabels({'0','1','10','100'});
            ylabel('Phosphorylation(AU)','FontSize',14);
            xlabel('time','FontSize',14);
            title(tit(i),'FontSize',14);
        end
end