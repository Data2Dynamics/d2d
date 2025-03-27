close all;
arInit;
global ar;

model = 'combined_model';
arSetTitle( ['Compiling ', model, ' ...'] );

addpath('/home/jvanlier/D2D/d2d/arFramework3/');
arInit;

arLoadModel(model);

loadData = 1;
if ( loadData )
    cull = {'time', @(t)(str2double(t)>240), 'input_dorso', @(dorso)(str2double(dorso)>0)};
    
    % qPCR data
    arLoadData( 'qpcr/qpcr_all',                    1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'qpcr/qpcr_APAP',                   1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'qpcr/qpcr_DR_BMP',                 1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'qpcr/qpcr_DR_IL6',                 1, 'csv', true, 'RemoveConditions', cull );     % OK

    % Western blot data
    arLoadData( 'wb_data/STAT3_TC',                 1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'wb_data/STAT3_TC_APAP',            1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'wb_data/SMAD_DR',                  1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'wb_data/STAT3_DR_IL6',             1, 'csv', true, 'RemoveConditions', cull );     % OK
    
    arLoadData( 'gp130/gp130_degradation',          1, 'csv', true );                               % OK
    arLoadData( 'SOCS3_wb/Anja_2018_13_SOCS3',      1, 'csv', true, 'RemoveConditions', cull );     % OK
    arLoadData( 'SOCS3_wb/Anja_2019_03_SOCS3',      1, 'csv', true, 'RemoveConditions', cull );     % OK
    
    arLoadData( 'BSA/protein_conc',                 1, 'csv', true );                               % OK
    
    arLoadData( 'qpcr/qpcr_decay',      1, 'csv', true, 'RemoveConditions', {'time', @(t)(str2double(t)<60||(str2double(t)==150))} );
    
    arLoadData('mlc/mlc',1,'xls',true);
end

%% Mock data file used to equilibrate the system
arLoadData('steadystate/steadystate_dcf_w_smad',1,'csv',true)               % OK

%% Compile
arCompileAll;

%% Model requires tight equilibration tolerances
ar.config.eq_tol = 1e-10;
ar.config.rtol = 1e-10;
ar.config.atol = 1e-10; 
edit Set
%% Pre-equilibration settings (use the steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;

%% Optimization settings
ar.config.optimizer = 5;
ar.config.optim.TolX = 1e-8;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.maxsteps = 100*ar.config.maxsteps;

%% Pre-equilibration settings (use the steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

arDisableData(arFindData('steadystate','names'));

% Remove the extremely high dose of APAP since these cells were of questionable viability
removeAPAP50;

% Remove a few datapoints that were identified as outliers
removeBad;

fixSelectedParameters;

arLink;

save(model, 'ar');