% arConditon2NewDataStruct(m,c,x,makeParamsSpecific,mtarget, useRealisticTimes)
% 
%   This function creates a new data struct from info given in a specific
%   condition. It was established to fit a function to the ODE solutions.
% 
% 
%   m   model index
%   c   condition index
%   x   x index or x name 
%   mtarget
%   useRealisticTimes   [true] CalculateRealisticTimeGrid is called
%                       false: tFine statt CalculateRealisticTimeGrid
% 
% 
%   ar.model(m).condition(c).xExpSimu(x) is set as data.yExp
% 
%   In this implementation, the transient function has to be added to have
%   ar.model.fy, ar.model.fystd etc 
% 
% Example:
% arConditon2NewDataStruct2(1,1,1)


function D = arConditon2NewDataStruct2(m,c,x,makeParamsSpecific,mtarget,useRealisticTimes)
if ~exist('useRealisticTimes','var') || isempty(useRealisticTimes)
    useRealisticTimes = true;
end
if ~exist('mtarget','var') || isempty(mtarget)
    mtarget = [];
end

global ar
if isnumeric(x)
    ix = x;
elseif ischar(x)
    ix = strmatch(x,ar.model(m).x,'exact');
end
if length(ix)~=1
    error('length(ix)~=1')
else
    x = ar.model(m).x{ix}; % be sure that it is a x name (char)
end

if ~exist('makeParamsSpecific','var') || isempty(makeParamsSpecific)
    makeParamsSpecific = true;
end

name = sprintf('Model%i_Cond%i_%s',m,c,x);
fprintf('Creating data struct %s ...\n',name);

if isempty(mtarget)
    mtarget = strmatch('TransientFunction_ForConditionFit2',{ar.model.name},'exact');
    if length(mtarget)~=1
        disp('In this implementation, the transient function has to be added to have ar.model.fy, ar.model.fystd etc ');
        error('Model ''TransientFunction_ForConditionFit2'' not available. Please create a model structure and add ''TransientFunction'' to it via arLoadModel.');
    end
end

tlim = [min(ar.model(m).condition(c).tFine),max(ar.model(m).condition(c).tFine)];
if useRealisticTimes
    tgrid = CalculateRealisticTimeGrid(ar.model(m).condition(c).tExp,tlim);
    xgrid = interp1(ar.model(m).condition(c).tFine,ar.model(m).condition(c).xFineSimu(:,ix),tgrid);
else
    tgrid = ar.model(m).condition(c).tFine;
    xgrid = ar.model(m).condition(c).xFineSimu(:,ix);
end

args = cell(0);
args{end+1} = 'tExp';       args{end+1} = tgrid;
args{end+1} = 'yExp';       args{end+1} = xgrid;
args{end+1} = 'yExpStd';    args{end+1} = NaN(size(xgrid));

% Not required, I'll prefer keeping it sparse:
% args{end+1} = 'yExpRaw';    args{end+1} = ar.model(m).condition(c).xExpSimu;
% args{end+1} = 'yExpStdRaw'; args{end+1} = NaN(size(ar.model(m).condition(c).xExpSimu));

args{end+1} = 'yNames';          args{end+1} = {x};
args{end+1} = 'y';          args{end+1} = {'TransientFunction2'};

% Not required: already specified in the model
% args{end+1} = 'fy';         args{end+1} = {'Signal_TransientFunction'}; 
% args{end+1} = 'fystd';      args{end+1} = {'sd_TF'};
args{end+1} = 'doseresponse'; args{end+1} = 0;

if rem(length(args),2)~=0
    error('arguments args has to be provided in pairs.')
end

% make parameters of the transient fucntion name-specific:
pold = {};
fp = {};
if makeParamsSpecific
    pold{end+1} = 'timescale_sust'; fp{end+1} = sprintf('timescale_sust_%s',name);
    pold{end+1} = 'timescale_trans'; fp{end+1} = sprintf('timescale_trans_%s',name);
    pold{end+1} = 'amp_sust'; fp{end+1} = sprintf('amp_sust_%s',name);
    pold{end+1} = 'amp_trans'; fp{end+1} = sprintf('amp_trans_%s',name);
    pold{end+1} = 'signum_TF'; fp{end+1} = sprintf('signum_TF_%s',name);
    pold{end+1} = 'offset_TF'; fp{end+1} = sprintf('offset_TF_%s',name);
    pold{end+1} = 'toffset_TF'; fp{end+1} = sprintf('toffset_TF_%s',name);
    pold{end+1} = 'sd_TF';     fp{end+1} = sprintf('sd_TF_%s',name);
end
pold{end+1} = 'maxt_TF'; fp{end+1} = ['(',num2str(range(ar.model(m).condition(c).tExp)),')'];
pold{end+1} = 'init____dummy___'; fp{end+1} = ['(0)'];

D = arCreateDataStruct(mtarget,pold,fp,args{:});



