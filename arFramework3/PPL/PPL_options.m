% options = PPL_options(varargin)
% 
% Set options for calculation of prediction profile likelihood. This function 
% must be called before use of arPPL. 
%
% The options are stored in ar.ppl.options.
% 
%   varargin        pairs optionname, optionvalue
%                   The option names are case INsensitive. 
% 
% The following options are available:
% 
%   Integrate       Integrate prediction/validation confidence line over
%                   time, or just compute profiles for specified t 
%                   [false]
%
%   doPPL           integrate prediction bands with assumed error on measurements
%                   [false]
%                   This default means that validation confidence intervals
%                   are calculated instead of prediction confidence
%                   intervals. The former gives CIs wrt data with non-zero
%                   standard deviation.
%
%   n_steps_profile  total steps in computation of full profile per side 
%                    [50]
%
%   xstd            standard deviation of additional data point. The used
%                   algorithm will try to set this automatically in arPPL,
%                   but if this is not possible it will be set to
%                   [0.1]
%                   Note that if prediction profiles are to be calculated
%                   the choice of xstd is crucial. Values smaller than the
%                   width of the prediction profile are recommended (which
%                   can only be found after the calculation has been
%                   attempted once before).
%
%   gamma           correction strength for computation of integration
%                   [1]
%
%   stepsize        Integration steps in time, e.g. 1, 0.5 ... 
%                   [default is tFine stepsize]
%
%   tEnd            optional stopping time for the confidence band integration
%                   [default is tLim(end)]
%
%   rel_increase    increase of xTrial in computation of full profile for t's
%                   [default is 0.15]
%
%   ed_steps        take integration steps or go parallel to orginial state
%                   solution 
%                   [true] 
%
%   fineInt         do precise integration with a correction cycle at every
%                   integration step 
%                   [true]
%
%   backward        perform integration backwards in time 
%                   [false]
%
%   whichT          which entry of the provided time-vector is taken as
%                   starting point for integration (argument is the
%                   corresponding index) [1]
%
%   dir             1 for only upper CI profile
%                   -1 for lower, 0 for both 
%                   [default is 0]
%                   Doesn't work in the case of ('integrate',false).
%
%   alpha_level     Confidence level for profile 
%                   [default 0.05]
%
% Helge Hass, 2014 (helge.hass@fdm.uni-freiburg.de)
% Written by Helge, tried to be documented by Clemens/Tim.
% 
% See also: arPPL
% 
% Example:
% >> PPL_options('alpha_level',.02,'tend',10)
% >> ar.ppl.options
% >> arPPL

function options = PPL_options(varargin)
global ar
if (nargin == 0) && (nargout == 0)
  fprintf(['     Integrate: [ Integrate prediction profiles (true) or ' ...
           '| no integration: {false} ]\n']);
  fprintf('          Gamma_factor: [ correction strength {1} ]\n');
  fprintf('        stepsize: [ integration step size {tFine} ]\n');  
  fprintf('            tEnd: [ Alternative integration endpoint ]\n'); 
  fprintf('    rel_increase: [ %% of x increase in profile calculation {0.15} ]\n');
  fprintf(['        ed_steps: [ Integration {true} | or parallel to initial trajectory: ' ...
           'false ]\n']);
  fprintf('         fineInt: [ Correction cycle after each step {true} | false ]\n');  
  fprintf('        backward: [ Integrate backwards in time: true | {false} ]\n');
  fprintf('            xstd: [ standard dev of auxiliary data point {auto calculation} ]\n');
  fprintf('           doPPL: [ true = prediction or false = validation profiles {false} ]\n');
  fprintf('         n_steps_profile: [ steps of profile in each direction {50} ]\n');
  fprintf('          whichT: [ seed for integration (in time vector) {1} ]\n');
  fprintf('             dir: [ only upper(1)/lower(-1) profile or both {0}? ]\n');
  fprintf('     alpha_level: [ Confidence level {0.05} ]\n');
  fprintf('\n');
  return;
end

% all possible option names:
Names = [  
    'gamma           '
    'stepsize        '
    'integrate       '
    'tEnd            '
    'rel_increase    '
    'ed_steps        '             
    'fineInt         '
    'backward        '
    'xstd            '
    'doPPL           '
    'n_steps_profile ' 
    'whichT          '
    'dir             '
    'alpha_level     '
    ];

%Set some default options
if(~isfield(ar.ppl,'options'))
    ar.ppl.options.integrate = false;
    ar.ppl.options.rel_increase = 0.15;
    ar.ppl.options.ed_steps = true;
    ar.ppl.options.fineInt = true;
    ar.ppl.options.backward = false;
    ar.ppl.options.xstd = 0.1;
    ar.ppl.options.doPPL = false;
    ar.ppl.options.n_steps_profile = 50;
    ar.ppl.options.whichT = 1;
    ar.ppl.options.dir = 0;
    ar.ppl.options.alpha_level = 0.05;
    ar.ppl.options.gamma = 1;
end

m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
% options = [];
% for j = 1:m
%   options.(deblank(Names(j,:))) = [];
% end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(message('No proper name %s', arg));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        ar.ppl.options.(deblank(Names(j,:))) = val;
        if(strcmp(deblank(Names(j,:)),'xstd'))
            ar.ppl.xstd_auto = 0;
        end
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error(message('No argument for a given option'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error(message('Option must be a char', arg));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error('Could not find option %s', arg);
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
            matches = deblank(Names(j(1),:));
        for k = j(2:length(j))'
                matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
        end
            error(message('Found more than one match for option %s',arg));
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    ar.ppl.options.(deblank(Names(j,:))) = arg;
    if(strcmp(deblank(Names(j,:)),'xstd'))
        ar.ppl.xstd_auto = 0;
    end
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(message('No value found for option %s', arg));
end

options = 'Did set PPL options!';
