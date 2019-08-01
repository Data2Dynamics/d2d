% This function approximates the dynamics of all model conditions and
% states.

function fits = arApproximateTimeCoursesByTransientFunction
global ar

if isempty(which('arFitTransient'))
    error('Please add folder to the library for the transient function, \ne.g. via %s','addpath(''...\d2d\arFramework3\Examples\ToyModels\TransientFunction\TransientFunction_library'')');
end

%% Step 1: Create data structs from conditions:
def_file = [fileparts(which('arInit')),filesep,'Examples',filesep,'ToyModels',filesep,'TransientFunction',filesep,'TransientFunction_library',filesep,'TransientFunction_ForConditionFit.def'];
copyfile(def_file,'Models');
arLoadModel('TransientFunction_ForConditionFit');
fits = cell(0);
for m=1:length(ar.model)
    for ix=1:length(ar.model(m).x)
        for c=1:length(ar.model(m).condition)
            fits{end+1} = struct;
            % since each condition is fitted individually below, no condidtoin-specific prameters are required.
            fits{end}.data = arConditon2NewDataStruct(m,c,ix,false);
            fits{end}.u = ar.model(end).u;
            fits{end}.m = m;
            fits{end}.c = c;
            fits{end}.ix = ix;
            fits{end}.x = ar.model(m).x{ix};
            fits{end}.condition = struct;
            fits{end}.condition.pold = ar.model(m).condition(c).pold;
            fits{end}.condition.fp = ar.model(m).condition(c).fp;
            fits{end}.condition.p = ar.model(m).condition(c).p;
            fits{end}.label = sprintf('Model%i_%s_Condition%i',fits{end}.m,fits{end}.x,fits{end}.c);
        end
    end
end
fprintf('%i data-structs created. Start fitting now...\n',length(fits));
ar.model = ar.model(1:(end-1)); % remove transient function

%% Step2: Fit the transient functions to each data struct:

for d=1:length(fits)
    arInit;
    arLoadModel('TransientFunction_ForConditionFit');
    arAddDataStruct(fits{d}.data);
    arCompileAll;
    
    boundfactor = 2;% In this setting (at least 20 data points, the following bounds are better)
    Initialize_FitTransient(boundfactor);
    ar.p(ar.qFit==1) = (ar.fit_transient.bounds.ub(ar.qFit==1)+ar.fit_transient.bounds.lb(ar.qFit==1))/2;
    arFitTransient
    %     arFitLHS(5)  % 5 fits should be enough for 7 parameters!
    
    arSimu(false,true,true);
    
    fits{d}.p = ar.p;
    fits{d}.merit = arGetMerit;
    fits{d}.pLabel = ar.pLabel;
    fits{d}.tFine = ar.model.data.tFine;
    fits{d}.yFineSimu = ar.model.data.yFineSimu;
    fits{d}.ystdFineSimu = ar.model.data.ystdFineSimu;
    
    fprintf('Fit %4i out of %5i done.\n',d,length(fits));
    if rem(d,10)==0
        save fits_tmp fits
    end
    
%     arPlot
%     drawnow
%     pause(0.1)
%     close all
end
