function opts = makeOptions(int, varargin)
% Create a JEInterface options set from scratch. Possible fields are:
%'Display';
%'MaxFunEvals';'MaxIter';'TolFun';'TolFunEvals';TolX';'TolXEvals', which,
% except for 'TolFunEvals' and 'TolXEvals', are used similar to the optimset. 
% The latter two are interpreted as the numbers of evaluations required
% to assume convergence. Default values are TolXEvals=TolFunEvals=200,
% TolX=TolFun=1e-4, MaxFunEvals uses a default from EvA2.
% Further options are 'CreateStopBox' and 'NiceSleepTime' which are
% JE-specific. The 'CreateStopBox' option toggles the GUI box with stop
% button provided for graceful interruption of an optimization run.
% In cases without X-forwarding, a deactivation by setting it to 0 may be required,
% but note that killing optimization may break the threading mechanism. Try stopOptimize with
% the 'kill' option.
%
% 'NiceSleepTime' allows setting a sleep time (in ms) for the concurrent
% threading. Technically, either the Matlab
% or the Java thread is always waiting for the other one, which may require 
% full CPU if no sleep time is employed. However, for quick function evaluations,
% even minor sleep times substantially reduce allover efficiency, so we advice
% to set small sleep times (e.g. 5 ms) only if a single function evaluation 
% takes considerably longer and leave the sleep time at zero otherwise.
%
% Notice that this method creates a parameter set but does not assign it
% to the interface instance. Use setOptions to do that.

allfields = {'Display'; 'MaxFunEvals';'MaxIter';'TolFun';'TolFunEvals';...
    'TolX'; 'TolXEvals'; 'CreateStopBox'; 'NiceSleepTime'};

nvararg=nargin-1;
if rem(nvararg,2)==1
    error('Pass options in name-value pairs!');
end

% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
opts = struct(structinput{:});

stdSet=optimset();

% standard options:
opts.('MaxFunEvals') = eva2.OptimizerFactory.getDefaultFitCalls;
opts.('TolX') = 1e-4;
opts.('TolXEvals') = 200;
opts.('TolFun') = 1e-4;
opts.('TolFunEvals') = 200;
opts.('CreateStopBox')=1;
opts.('NiceSleepTime')=0;

for i=1:nvararg/2
    name=varargin{2*i-1};
    value=varargin{(2*i)};
    % parse arguments
    if ~ischar(name)
        error('Expected char parameter name at index %d!', 2*i+1);
    else
        optIndex=strmatch(name,allfields, 'exact');
        if isempty(optIndex)
            error('Unknown option %s !', name);
        else
            switch name
                case {'TolFunEvals', 'TolXEvals'} % for special fields, integers > 1 are allowed
                    if ~isscalar(value) || ~isnumeric(value) || round(value)<1
                        error('Invalid value type for %s, expecting positive numeric scalar!', name);
                    end;
                    value=round(value);
                case 'CreateStopBox'
                    if ~isscalar(value) || ~isnumeric(value) || round(value)<0 || round(value)>1
                        error('Invalid value type for %s, expecting 0 or 1!', name);
                    end;
                    value=round(value);
                case 'NiceSleepTime'
                    if ~isscalar(value) || ~isnumeric(value) || round(value)<0 ;
                        error('Invalid value type for %s, expecting numeric scalar >= 0!', name);
                    end;
                    value=round(value);
                otherwise                % test using optimset
                    optimset(stdSet, name, value);
            end
            % assign to struct
            opts.(allfields{optIndex,:}) = value;
        end
    end
end
