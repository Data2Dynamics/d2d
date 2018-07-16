% delta = arCompareConfig
% 
%   config_field    which configuration of ar.config should be altered.
%                   It works also if serveral struct levels are specified.
% 
%                   Examples: 'useEvents', 'optimizer', 'optim.Maxiter'...
% 
%    config_vals    Either specified as numeric array or as cells
% 
%    userFun        An optional user-defined function which is applied
%                   after configs are altered. This function can change
%                   global ar.
% 
% This function evaluates the impace of ar.config.useEvents at the level of
% residuals, xExpSimu and yExpSimu
% 
% 
% Example:
% delta = arCompareConfig('useEvents',[0,1])
% delta = arCompareConfig('useEvents',[0,1],@arFit)
% delta = arCompareConfig('rtol',logspace(-8,-5,4),@arFit);
% 

function delta = arCompareConfig(config_field,config_vals,user_fun)

global ar

if ~exist('config_field','var') || isempty(config_field)
   error('Please provide a config field, e.g. ''useEvents''')
elseif ~ischar(config_field)
    error('~ischar(config_field)');
end
if ~exist('config_vals','var') || isempty(config_vals)
   error('Please provide a config values, e.g. [0,1]')
elseif isnumeric(config_vals)
    config_vals = num2cell(config_vals);
elseif ~iscell(config_vals)
    error('Please specify config_vals as numerical array or as cell array.')
end
if ~exist('user_fun','var') || isempty(user_fun)
    user_fun = '';
end


% remember old value
ar0 = arDeepCopy(ar);

try
    ars = cell(size(config_vals));
    for ii=1:length(config_vals)
        ar = arDeepCopy(ar0); % reset original state
        evstr = sprintf('ar.config.%s = config_vals{ii};',config_field);
        eval(evstr)
        
        arCheckCache(true); % invalidate cache
        arCalcMerit
        arSimu(false,false,true)
        if ~isempty(user_fun)
            feval(user_fun);
        end
        ars{ii} = arDeepCopy(ar);
        ars{ii}.merit = arGetMerit(true);
    end
    
    % compare
    delta.maxabs = struct;
    delta.maxabs.merit = -Inf;
    delta.maxabs.res = -Inf;
    delta.maxabs.xExpSimu = -Inf;
    delta.maxabs.yExpSimu = -Inf;
    
    for i1=1:(length(ars)-1)
        for i2=(i1+1):length(ars)
            fn = sprintf('compare_%i_vs_%i',i1,i2);
            delta.(fn) = struct;

            delta.(fn).merit = ars{i2}.merit-ars{i1}.merit;
            delta.maxabs.merit = max(delta.maxabs.merit,max(abs(delta.(fn).merit)));            
            
            delta.(fn).res = ars{i2}.res-ars{i1}.res;
            delta.maxabs.res = max(delta.maxabs.res,max(abs(delta.(fn).res)));
            
            delta.(fn).xExpSimu = cell(size(ars{i2}.model));
            for m=1:length(ars{i2}.model)
                delta.(fn).xExpSimu{m} = cell(size(ars{i2}.model(m).condition));
                for c=1:length(ars{i2}.model(m).condition)
                    delta.(fn).xExpSimu{m}{c} = ars{i2}.model(m).condition(c).xExpSimu - ars{i1}.model(m).condition(c).xExpSimu;
                    if ~isempty(delta.(fn).xExpSimu{m}{c})
                        delta.maxabs.xExpSimu = max(delta.maxabs.xExpSimu,max(abs(delta.(fn).xExpSimu{m}{c}(:))));
                    end
                end
            end
            
            
            delta.(fn).yExpSimu = cell(size(ars{i2}.model));
            for m=1:length(ars{i2}.model)
                delta.(fn).yExpSimu{m} = cell(size(ars{i2}.model(m).data));
                for d=1:length(ars{i2}.model(m).data)
                    delta.(fn).yExpSimu{m}{d} = ars{i2}.model(m).data(d).yExpSimu - ars{i1}.model(m).data(d).yExpSimu;
                    if ~isempty(delta.(fn).yExpSimu{m}{c})
                        delta.maxabs.yExpSimu = max(delta.maxabs.yExpSimu,max(abs(delta.(fn).yExpSimu{m}{d}(:))));
                    end
                end
            end
        end
    end
    
catch ERR
    ar = arDeepCopy(ar0); % reset original state

    rethrow(ERR)
end

ar = arDeepCopy(ar0); % reset original state

disp('----------Result:-----------')
fprintf('\nMax. absolute differences (might depend on ar.p):\n')
fprintf('  xExpSimu:  %f\n',delta.maxabs.xExpSimu);
fprintf('  yExpSimu:  %f\n',delta.maxabs.yExpSimu);
fprintf('  residuals: %f\n',delta.maxabs.res);
fprintf('  merit-fkt: %f\n',delta.maxabs.merit);

if delta.maxabs.res>0.1
    disp('strong impact => be very careful !!')
elseif delta.maxabs.res>1e-3
    disp('relevant impact => be careful!')
elseif delta.maxabs.res<1e-6
    disp('rather weak impact for the current parameter set.')
end
disp('----------------------------')





        

