% delta = arCompareUseEvents
% 
% This function evaluates the impace of ar.config.useEvents at the level of
% residuals, xExpSimu and yExpSimu
%
%   delta = 
%                 res: [1x1428 double]
%          res_maxabs: 14.6769
%     xExpSimu_maxabs: 0.7119
%            xExpSimu: {{1x109 cell}}
%     yExpSimu_maxabs: 1.2539e+03
%            yExpSimu: {{1x123 cell}}    
% 
% Example:
% delta = arCompareUseEvents
% 
% 94 event times found over 109 conditions.
% -2*log(L) = 18511.4, 713 data points, 46 free parameters, , 4.29422 violation of 26 prior assumptionsdata chi^2 = 11203.1
% -2*log(L) = 7995.25, 713 data points, 46 free parameters, , 4.29422 violation of 26 prior assumptionsdata chi^2 = 687
% Absolute differences between useEvents=0 vs. useEvents=1 (only for current ar.p):
%   xExpSimu:  14.676856
%   yExpSimu:  14.676856
%   residuals: 14.676856
% useEvents flag has a strong impact => useEvents=1 strongly suggested.
%        

function delta = arCompareUseEvents

global ar

% first check, whether events are present:
tev = [];
ncond = 0;
for m=1:length(ar.model)
    ncond = ncond + length(ar.model(m).condition);
    for c=1:length(ar.model(m).condition)
        if isfield(ar.model(m).condition(c),'tEvents')
            tev = [tev,setdiff(ar.model(m).condition(c).tEvents,ar.model(m).condition(c).tstart)];
        end
    end
end
tev

disp('----------------------------')

if isempty(tev)
    disp('No events found. Call arFindInputs or arAddEvents first.')
    delta = [];
    return
else
    fprintf('%i event times found over %i conditions.\n',length(tev),ncond);
end

useEv0 = ar.config.useEvents;

try
    ar.config.useEvents = 0;
    disp('ar.config.useEvents=0:')
    arCheckCache(true); % invalidate cache
    arCalcMerit
    arSimu(false,false,true)
    ar0 = arDeepCopy(ar);
    
    disp('ar.config.useEvents=1:')
    ar.config.useEvents = 1;
    arCheckCache(true); % invalidate cache
    arCalcMerit
    arSimu(false,false,true)
    
    % compare
    delta.res = ar.res-ar0.res;
    delta.res_maxabs = max(abs(delta.res));
    
    delta.xExpSimu_maxabs = -Inf;
    delta.xExpSimu = cell(size(ar.model));
    for m=1:length(ar.model)
        delta.xExpSimu{m} = cell(size(ar.model(m).condition));
        for c=1:length(ar.model(m).condition)
            delta.xExpSimu{m}{c} = ar.model(m).condition(c).xExpSimu - ar0.model(m).condition(c).xExpSimu;
            if ~isempty(delta.xExpSimu{m}{c})
                delta.xExpSimu_maxabs = max(delta.xExpSimu_maxabs,max(abs(delta.xExpSimu{m}{c}(:))));
            end
        end
    end
    
    
    delta.yExpSimu_maxabs = -Inf;
    delta.yExpSimu = cell(size(ar.model));
    for m=1:length(ar.model)
        delta.yExpSimu{m} = cell(size(ar.model(m).data));
        for d=1:length(ar.model(m).data)
            delta.yExpSimu{m}{d} = ar.model(m).data(d).yExpSimu - ar0.model(m).data(d).yExpSimu;
            if ~isempty(delta.yExpSimu{m}{c})
                delta.yExpSimu_maxabs = max(delta.yExpSimu_maxabs,max(abs(delta.yExpSimu{m}{d}(:))));
            end
        end
    end
    
catch ERR
    ar.config.useEvents = useEv0;
    rethrow(ERR)
end

ar.config.useEvents = useEv0;

fprintf('\nMax. absolute differences between useEvents=0 vs. useEvents=1 (only for current ar.p):\n')
fprintf('  xExpSimu:  %f\n',delta.xExpSimu_maxabs);
fprintf('  yExpSimu:  %f\n',delta.yExpSimu_maxabs);
fprintf('  residuals: %f\n',delta.res_maxabs);

if delta.res_maxabs>0.1
    disp('useEvents flag has a strong impact => useEvents=1 strongly suggested.')
elseif delta.res_maxabs>1e-3
    disp('useEvents flag has a relevant impact => useEvents=1 suggested.')
elseif delta.res_maxabs<1e-6
    disp('useEvents flag has a weak impact for the current parameter set.')
end
disp('----------------------------')





        

