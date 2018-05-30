datafiles% Load data set to next free slot
%
% arLoadData(name, m, extension, removeEmptyObs, opts)
%
% name                      filename of data definition file
% m                         target position (int) for model                [last loaded model]
%                           or: model name (string)
% extension                 data file name-extension: 'xls', 'csv'         ['xls']
%                           'none' = don't load data                           
% removeEmptyObs            remove observation without data                [false]
% opts                      additional option flags
%
% optional option flags are:
% 'RemoveConditions'        This flag followed by a list of conditions will
%                           allow you to filter the data that you load. Note
%                           that the function takes both values, strings
%                           or function handles. When you provide a function
%                           handle, the function will be evaluated for each
%                           condition value (in this case input_dcf). You
%                           should make the function return 1 if the condition
%                           is to be removed. Note that the input to the
%                           function is a *string* not a number.
%                           Example:
%                            arLoadData( 'mydata', 1, 'csv', true, 'RemoveConditions', ...
%                            {'input_il6', '0', 'input_dcf', @(dcf)str2num(dcf)>0};
% 'RemoveEmptyConds'        This flag allows you to remove conditions that have
%                           no data points.
% 'ResampleDoseResponse'    This flag allows you to increase the resolution
%                           of dose responses for plotting purposes. Note
%                           that this makes model evaluation and
%                           compilation much more time intensive (hence it
%                           is off by default).
% 'ResamplingResolution'    Number of extra points to interpolate dose
%                           response (default: 25)
% 'RefineLog'               Perform the refinement on a log scale
% 'expsplit'                Split replicate specific parameters based on the 
%                           column identifier specified in the next argument
%                           Example:
%                               arLoadData('mRNA/mRNA_pretreatment', 1, 'csv', true, 'expsplit', 'nExpID');
%                           will recognize parameters with nExpID in the name 
%                           and replace them with nExpID0, nExpID1 etc. if
%                           the excel data file has a column named nExpID  
%                           which specifies these values.
%
% 'RemoveObservables'   This can be used to omit observables. Simply pass a
%                       cell array with names of observables that should be 
%                       ignored in the data.
%
% 'IgnoreInputs'        Ignore input overrides. Simply pass a cell array
%                       with names of inputs that should be ignored (the 
%                       model default will be used instead).
% 
% 'DataPath'            Path to the data files.
%                       Default: DataPath = 'Data/'
% 
%
% The data file specification is as follows:
%
%   In the first column
%   1)    Measurement time points (are allowed to occur multiple times).
%
%   In the following columns (in any order):
%   2)    Experimental conditions (e.g. "input_IL6" and "input_IL1").
%   3)    The data points for the individual observables (e.g. "P_p38_rel").
%
%   Note:
%   1)    No mathematical symbols are allowed in the column headers (e.g. "+")
%   2)    I have always used input_ as a prefix for stimulations. Regarding
%         observables, the suffixes "_rel" and "_au" refer to relative
%         phosphorylation and arbitrary units.
%
% Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)

function arLoadData(name, m, extension, removeEmptyObs, varargin)

global ar

arFprintf( 3, 'Parsing input arguments...\n' );

if(isempty(ar))
    error('please initialize by arInit')
end

switches = { 'dppershoot', 'removeconditions', 'removeobservables', 'splitconditions',...
    'removeemptyconds', 'expsplit', 'resampledoseresponse', 'resamplingresolution',...
    'refinelog', 'ignoreinputs', 'detectionlimit', 'datapath'};
extraArgs = [ 1, 1, 1, 1, ...
    0, 1, 0, 1, ...
    0, 1, 1, 1 ];
description = { ...
    {'', 'Multiple shooting on'} ...
    {'', 'Ignoring specific conditions'} ...
    {'', 'Ignoring specific observables'} ...
    {'', 'Split data set into specific conditions'}, ...
    {'', 'Removing conditions without data'}, ...
    {'', 'Splitting conditions by specific data column'}, ...
    {'', 'Resampling dose response'}, ...
    {'', 'Resampling with custom resolution'}, ...
    {'', 'Resampling on log scale'}, ...
    {'', 'Ignoring specific inputs'}, ...
    {'', 'Working with detection limit'}, ...
    {'', 'Path to the data files'}};
    
opts = argSwitch( switches, extraArgs, description, 1, varargin );
if isempty(opts.datapath_args)
    DataPath = 'Data/';
else
    DataPath = opts.datapath_args;
    if DataPath(end)~='/' && DataPath(end)~='\'
        DataPath = [DataPath,'/'];
    end
end


% load model from mat-file
if(~exist(DataPath,'dir'))
    error('folder %s does not exist',DataPath)
end
if strcmp(strrep(name,' ',''),name)~=1
    name
    error('File names should not contain empty spaces. Please remove it.');
end
if(~exist([DataPath, name, '.def'],'file'))
    if(~exist([DataPath, name '.xls'],'file') && ~exist([DataPath, name '.csv'],'file') && ~exist([DataPath, name '.xlsx'],'file'))
        error('data definition file %s.def does not exist in folder %s', name,DataPath)
    else
        arFprintf(1, '\ncreating generic .def file for %s/%s ...\n', DataPath, name);
        copyfile(which('data_template.def'),[DataPath, name, '.def']);
    end
else
    if(~exist([DataPath, name '.xls'],'file') && ~exist([DataPath, name '.csv'],'file') && ~exist([DataPath, name '.xlsx'],'file'))
        warning('data file corresponding to %s.def does not exist in folder %s/', DataPath, name)
    end
end
    
if(~exist('m','var') || isempty(m))
    m = length(ar.model);
end
if(exist('m','var') && ischar(m))
    for jm=1:length(ar.model)
        if(strcmp(m, ar.model(jm).name))
            m = jm;
        end
    end
    if(ischar(m))
        error('Model %s was not found', m);
    end
end

if(exist('extension','var') && isnumeric(extension) && ...
        ~(isempty(extension) && nargin>3))
    error(['arLoadData(name, m, d, ...) input argument d is deprecated !!! ' ...
        'Please see new usage arLoadModel(name, m, extension, removeEmptyObs) and function help text.']);
end

if(isfield(ar.model(m), 'data'))
    d = length(ar.model(m).data) + 1;
else
    ar.model(m).data = [];
    d = 1;
end

if(~exist('extension','var') || isempty(extension))
    extension = 'xls';
    
    % auto-select extension if not specified
    if exist([DataPath, name '.xlsx'],'file')
        extension = 'xlsx';
    elseif exist([DataPath, name '.xls'],'file')
        extension = 'xls';
    elseif exist([DataPath, name '.csv'],'file')
        extension = 'csv';
    end
end
if(~exist('removeEmptyObs','var'))
    removeEmptyObs = false;
else
    if(ischar(removeEmptyObs))
        error(['arLoadData(name, m, d, ...) input argument d is deprecated !!! ' ...
            'Please see new usage arLoadModel(name, m, extension, removeEmptyObs) and function help text.']);
    end
end


%%
if ( opts.resampledoseresponse )
    if ( ~isnumeric( opts.resamplingresolution_args ) || isempty( opts.resamplingresolution_args ) )
        opts.resamplingresolution = 25;
    else
        opts.resamplingresolution = opts.resamplingresolution_args(1);
    end
end

if( opts.dppershoot )
    if( opts.dppershoot_args>0 )
        if(~isfield(ar,'ms_count_snips'))
            ar.model(m).ms_count = 0;
            ar.ms_count_snips = 0;
            ar.ms_strength = 0;
            ar.ms_threshold = 1e-5;
            ar.ms_violation = [];
        end
        dpPerShoot = opts.dppershoot_args;
    end
else
    dpPerShoot = 0;
end

% initial setup
ar.model(m).data(d).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
ar.model(m).data(d).path = [pwd,filesep,DataPath];

ar.model(m).data(d).uNames = {};

arFprintf(1, '\nloading data #%i, from file %s%s.def...', d, DataPath, name);

% Disable this if you are having problems because of the preprocessor
preprocessor = 1;
if ( ~preprocessor )
    fid = fopen([DataPath, name, '.def'], 'r');
else
    % Load into a struct
    fid.fn  = [DataPath name '.def'];
    fid.str = fileread([DataPath name '.def']);
    fid.pos = 1;
    arFprintf( 3, 'Running preprocessor...' );
    fid = arPreProcessor(fid);
    arFprintf( 3, ' [ OK ]\n' );
end

% DESCRIPTION
arFprintf( 3, 'Start parsing...' );
[str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(~strcmp(str{1},'DESCRIPTION'))
    arParsingError( fid, 'parsing data %s for DESCRIPTION', name);
end

% check version
if(strcmp(str{1},'DESCRIPTION'))
    % def_version = 1;
elseif(strcmp(str{1},'DESCRIPTION-V2'))
    arParsingError( fid, 'DESCRIPTION-V2 not supported yet');
else
    arParsingError( fid, 'invalid version identifier: %s', cell2mat(str{1}));
end

% read comments
arFprintf( 3, '[ OK ]\nReading description' );
[str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
ar.model(m).data(d).description = {};
while(~strcmp(str{1},'PREDICTOR') && ~strcmp(str{1},'PREDICTOR-DOSERESPONSE'))
    arFprintf( 3, '.' );
    ar.model(m).data(d).description(end+1,1) = str{1}; %#ok<*AGROW>
    [str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
end

% PREDICTOR
arFprintf( 3, '[ OK ]\nReading predictor...\n' );
if(strcmp(str{1},'PREDICTOR-DOSERESPONSE'))
    ar.model(m).data(d).doseresponse = true;
    [str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
    ar.model(m).data(d).response_parameter = cell2mat(str{1});
    arFprintf(2, 'dose-response to %s\n', ar.model(m).data(d).response_parameter);
else
    ar.model(m).data(d).doseresponse = false;
    ar.model(m).data(d).response_parameter = '';
    arFprintf(2, '\n');
end
[C, fid] = arTextScan(fid, '%s %s %q %q %n %n %n %n\n',1, 'CommentStyle', ar.config.comment_string);

ar.model(m).data(d).t = cell2mat(C{1});
ar.model(m).data(d).tUnits(1) = C{2};
ar.model(m).data(d).tUnits(2) = C{3};
ar.model(m).data(d).tUnits(3) = C{4};
ar.model(m).data(d).tLim = [checkNum(C{5}, 0) checkNum(C{6}, 10)];
ar.model(m).data(d).tLimExp = [checkNum(C{7}, ar.model(m).tLim(1)), checkNum(C{8}, ar.model(m).tLim(2))];

% INPUTS
arFprintf( 3, 'Reading inputs...' );
[str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(~strcmp(str{1},'INPUTS'))
    arParsingError( fid, 'parsing data %s for INPUTS', name);
end
[C, fid] = arTextScan(fid, '%s %q %q\n',1, 'CommentStyle', ar.config.comment_string);
ar.model(m).data(d).fu = ar.model(m).fu;
while(~strcmp(C{1},'OBSERVABLES'))
    arFprintf( 3, '.' );
    qu = ismember(ar.model(m).u, C{1}); %R2013a compatible
    if(sum(qu)~=1)
        arParsingError( fid, 'unknown input %s', cell2mat(C{1}));
    end
    
    % Ignore this replacement?
    ignoreInput = 0;
    if (opts.ignoreinputs)
        if ismember(C{1}, opts.ignoreinputs_args)
            ignoreInput = 1;
        end
    end
    if ( ~ignoreInput )
        % Input replacement description
        ar.model(m).data(d).fu(qu) = C{2};
        if(~isempty(cell2mat(C{3})))
            ar.model(m).data(d).uNames(end+1) = C{3};
        else
            ar.model(m).data(d).uNames{end+1} = '';
        end
    end
    [C, fid] = arTextScan(fid, '%s %q %q\n',1, 'CommentStyle', ar.config.comment_string);
end

% input parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fu, 'UniformOutput', false);
ar.model(m).data(d).pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}); %R2013a compatible

% OBSERVABLES
arFprintf( 3, '[ OK ]\nReading observables' );
if(isfield(ar.model(m),'y'))
    ar.model(m).data(d).y = ar.model(m).y;
    ar.model(m).data(d).yNames = ar.model(m).yNames;
    ar.model(m).data(d).yUnits = ar.model(m).yUnits;
    ar.model(m).data(d).normalize = ar.model(m).normalize;
    ar.model(m).data(d).logfitting = ar.model(m).logfitting;
    ar.model(m).data(d).logplotting = ar.model(m).logplotting;
    ar.model(m).data(d).fy = ar.model(m).fy;
else 
    ar.model(m).data(d).y = {};
    ar.model(m).data(d).yNames = {};
    ar.model(m).data(d).yUnits = {};
    ar.model(m).data(d).normalize = [];
    ar.model(m).data(d).logfitting = [];
    ar.model(m).data(d).logplotting = [];
    ar.model(m).data(d).fy = {};
end

[C, fid] = arTextScan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'ERRORS'))
    arFprintf( 3, '.' );
    qyindex = ismember(ar.model(m).data(d).y, C{1});
    if(sum(qyindex)==1)
        yindex = find(qyindex);
    elseif(sum(qyindex)==0)
        yindex = length(ar.model(m).data(d).y) + 1;
    else
        arParsingError( fid, 'multiple matches for %s', cell2mat(C{1}))
    end
    
    ar.model(m).data(d).y(yindex) = C{1};
    ar.model(m).data(d).yUnits(yindex,1) = C{2};
    ar.model(m).data(d).yUnits(yindex,2) = C{3};
    ar.model(m).data(d).yUnits(yindex,3) = C{4};
    ar.model(m).data(d).normalize(yindex) = C{5};
    ar.model(m).data(d).logfitting(yindex) = C{6};
    ar.model(m).data(d).logplotting(yindex) = C{6};
    ar.model(m).data(d).fy(yindex,1) = C{7};
    if(~isempty(cell2mat(C{8})))
        ar.model(m).data(d).yNames(yindex) = C{8};
    else
        ar.model(m).data(d).yNames(yindex) = ar.model(m).data(d).y(yindex);
    end
    [C, fid] = arTextScan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', ar.config.comment_string);
    if(sum(ismember(ar.model(m).x, ar.model(m).data(d).y{yindex}))>0) %R2013a compatible
        arParsingError( fid, '%s already defined in STATES', ar.model(m).data(d).y{yindex});
    end
    if(sum(ismember(ar.model(m).u, ar.model(m).data(d).y{end}))>0) %R2013a compatible
        arParsingError( fid, '%s already defined in INPUTS', ar.model(m).data(d).y{end});
    end
    if(sum(ismember(ar.model(m).z, ar.model(m).data(d).y{end}))>0) %R2013a compatible
        arParsingError( fid, '%s already defined in DERIVED', ar.model(m).data(d).y{end});
    end
    if(sum(ismember(ar.model(m).p, ar.model(m).data(d).y{end}))>0) %R2013a compatible
        arParsingError( fid, '%s already defined as parameter', ar.model(m).data(d).y{end});
    end
end

% observation parameters
arFprintf( 3, '[ OK ]\nComputing observation parameters...' );
varlist = cellfun(@symvar, ar.model(m).data(d).fy, 'UniformOutput', false);
ar.model(m).data(d).py = setdiff(setdiff(vertcat(varlist{:}), union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
if(isempty(ar.model(m).data(d).fy))
    arParsingError( fid, 'No OBSERVABLE specified. Specify an OBSERVABLE in the model or data definition file. See "Defining the OBSERVABLES".');
end
for j=1:length(ar.model(m).data(d).fy)
    varlist = symvar(ar.model(m).data(d).fy{j});
    ar.model(m).data(d).py_sep(j).pars = setdiff(setdiff(varlist, union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
    
    % exclude parameters form model definition
    ar.model(m).data(d).py_sep(j).pars = setdiff(ar.model(m).data(d).py_sep(j).pars, ar.model(m).px);
    ar.model(m).data(d).py_sep(j).pars = setdiff(ar.model(m).data(d).py_sep(j).pars, ar.model(m).pu);
end

% ERRORS
arFprintf( 3, '[ OK ]\nReading error models' );
if(isfield(ar.model(m),'y'))
    ar.model(m).data(d).fystd = ar.model(m).fystd;
else
    ar.model(m).data(d).fystd = cell(0);
end
[C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'INVARIANTS') && ~strcmp(C{1},'DERIVED') && ~strcmp(C{1},'CONDITIONS') && ~strcmp(C{1},'SUBSTITUTIONS'))
    arFprintf( 3, '.' );
    qyindex = ismember(ar.model(m).data(d).y, C{1});
    if ( qyindex == 0 )
        arParsingError( fid,  'Specified error model for non existent observable %s', C{1}{1} );
    end
    y_var_name = setdiff(symvar(ar.model(m).data(d).fy{qyindex}),ar.model(m).data(d).py);
    reg_string = ['((?<=\W)|^)(',C{1}{1},'|'];
    for jreg = 1:length(y_var_name)
        if(jreg<length(y_var_name))
            reg_string = [reg_string ,y_var_name{jreg},'|'];
        else
            reg_string = [reg_string ,y_var_name{jreg},')'];
        end
    end
    reg_string = [reg_string '((?=\W)|$)'];
    if(~isempty(regexp(C{2}{1},reg_string,'ONCE')) && ar.model(m).data(d).logfitting(qyindex))
       warning(['You are trying to set up a relative error model within a log transformation. \n%s' ...
        'Comment out this error if you want to proceed anyway. To implement an absolute error in log, \n' ...
        'you can try the approach: \nyObs = sd_yObs + 1/2 * (a+sqrt((a)^2)), a = (offset - yObs-sd_yObs) \n, with hard set or fitted offset (on log-scale) \n'],C{2}{1})
        arParsingError( fid, 'Revise error model')
    end
    if(sum(qyindex)==1)
        yindex = find(qyindex);
    elseif(sum(qyindex)==0)
        yindex = length(ar.model(m).data(d).y) + 1;
        warning('Specified error without specifying observation function (%s in %s). Proceed with caution!', C{1}{1}, ar.model(m).data(d).name);
    else
        arParsingError( fid, 'multiple matches for %s', cell2mat(C{1}))
    end
    ar.model(m).data(d).fystd(yindex) = C{2};
    [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
    
    if ( isempty( C{1} ) )
        arParsingError( fid, 'Unexpected end of file after error model (did you forget CONDITIONS?)')
    end
end

if(length(ar.model(m).data(d).fystd)<length(ar.model(m).data(d).fy))
    arParsingError( fid, 'some observables do not have an error model defined');
end

% Drop certain observables
if (opts.removeobservables)
    arFprintf( 3, '[ OK ]\nDropping specific observables...\n' );
    if ischar( opts.removeobservables_args )
        opts.removeobservables_args = {opts.removeobservables_args};
    end
    for a = 1 : length( opts.removeobservables_args )
        jo = 1;
        while( jo <= length( ar.model(m).data(d).y ) )
            jind = ismember( ar.model(m).data(d).y{jo}, opts.removeobservables_args );
            if ( sum(jind) > 0 )
                warning( '>> Explicitly removing %s!\n', ar.model(m).data(d).y{jo} );
                ar.model(m).data(d).y(jo) = [];
                ar.model(m).data(d).yUnits(jo,:) = [];
                ar.model(m).data(d).normalize(jo) = [];
                ar.model(m).data(d).logfitting(jo) = [];
                ar.model(m).data(d).logplotting(jo) = [];
                ar.model(m).data(d).fy(jo) = [];
                ar.model(m).data(d).yNames(jo) = [];
                ar.model(m).data(d).fystd(jo) = [];
            else
                jo = jo + 1;
            end
        end
    end
end

% error parameters
arFprintf( 3, 'Compute error parameters...' );
varlist = cellfun(@symvar, ar.model(m).data(d).fystd, 'UniformOutput', false);
ar.model(m).data(d).pystd = setdiff(vertcat(varlist{:}), union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ... %R2013a compatible
    ar.model(m).data(d).y), ar.model(m).t));
for j=1:length(ar.model(m).data(d).fystd)
    varlist = symvar(ar.model(m).data(d).fystd{j});
	ar.model(m).data(d).py_sep(j).pars = union(ar.model(m).data(d).py_sep(j).pars, ... %R2013a compatible
        setdiff(varlist, union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ar.model(m).data(d).y))); %R2013a compatible
    
    % exclude parameters form model definition
    ar.model(m).data(d).py_sep(j).pars = setdiff(ar.model(m).data(d).py_sep(j).pars, ar.model(m).px);
    ar.model(m).data(d).py_sep(j).pars = setdiff(ar.model(m).data(d).py_sep(j).pars, ar.model(m).pu);
end

% DERIVED
if(strcmp(C{1},'DERIVED'))
    arParsingError( fid, ['There is no need for a section DERIVED in data definition file! ' ...
        'Please remove and see usage in: ' ...
        'https://github.com/Data2Dynamics/d2d/wiki/Setting%20up%20models']);
end
% INVARIANTS
if(strcmp(C{1},'INVARIANTS'))
    arParsingError( fid, ['Section INVARIANTS in data definition file is deprecated! ' ...
        'Please remove and see usage in: ' ...
        'https://github.com/Data2Dynamics/d2d/wiki/Setting%20up%20models']);
end

% collect parameters needed for OBS
ptmp = union(ar.model(m).px, ar.model(m).pu);
ar.model(m).data(d).p = union(ptmp, union(ar.model(m).data(d).pu, ar.model(m).data(d).py)); %R2013a compatible
ar.model(m).data(d).pystd = setdiff(ar.model(m).data(d).pystd, ar.model(m).data(d).p); %Remove dynamic variables from error model parameters
ar.model(m).data(d).p = union(ar.model(m).data(d).p, ar.model(m).data(d).pystd); %R2013a compatible

% Union's behaviour is different when first arg is empty. In this case, a
% flip of the parameter vector is typically required.
if ( size( ar.model(m).data(d).p, 1 ) ~= 1 )
    ar.model(m).data(d).p = ar.model(m).data(d).p.';
end

% replace filename
ar.model(m).data(d).p = strrep(ar.model(m).data(d).p, '_filename', ['_' ar.model(m).data(d).name]);
ar.model(m).data(d).fy = strrep(ar.model(m).data(d).fy, '_filename', ['_' ar.model(m).data(d).name]);
ar.model(m).data(d).py = strrep(ar.model(m).data(d).py, '_filename', ['_' ar.model(m).data(d).name]);
ar.model(m).data(d).fystd = strrep(ar.model(m).data(d).fystd, '_filename', ['_' ar.model(m).data(d).name]);
ar.model(m).data(d).pystd = strrep(ar.model(m).data(d).pystd, '_filename', ['_' ar.model(m).data(d).name]);
for j=1:length(ar.model(m).data(d).py_sep)
    ar.model(m).data(d).py_sep(j).pars = strrep(ar.model(m).data(d).py_sep(j).pars, '_filename', ['_' ar.model(m).data(d).name]);
end

% SUBSTITUTIONS (beta)
substitutions = 0;
matVer = ver('MATLAB');
if ( strcmp(C{1},'SUBSTITUTIONS') )
    arFprintf( 3, '[ OK ]\nReading substitutions' );
    if(str2double(matVer.Version)>=8.4)
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
    else
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16);
    end    
    
    % Substitutions
    fromSubs = {};
    toSubs = {};
    ismodelpar = [];

    % Fetch desired substitutions
    while(~isempty(C{1}) && ~strcmp(C{1},'CONDITIONS'))
        arFprintf( 3, '.' );
        fromSubs(end+1)     = C{1}; %#OK<AGROW>
        toSubs(end+1)       = C{2}; %#OK<AGROW>
        ismodelpar(end+1)   = sum(ismember(ar.model(m).p, C{1})); %#OK<AGROW>

        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
    end

    if ( sum(ismodelpar) > 0 )
        s = sprintf( '%s\n', fromSubs{ismodelpar>0} );
        arParsingError( fid,  'Cannot substitute model parameters. These following parameters belong under CONDITIONS:\n%s', s );
    end

    % Perform selfsubstitutions
    arFprintf( 3, '[ OK ]\nPerforming self substitutions...' );
    if ( ~isempty(fromSubs) )
        substitutions = 1;
        toSubs = arSubsRepeated( toSubs, fromSubs, toSubs, str2double(matVer.Version) );
    end
    arFprintf( 3, '[ OK ]\n' );
end

% CONDITIONS
arFprintf( 3, 'Reading conditions' );
[C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);    
ar.model(m).data(d).fp = transpose(ar.model(m).data(d).p);
ptmp = ar.model(m).p;
qcondparamodel = ismember(ar.model(m).data(d).p, strrep(ptmp, '_filename', ['_' ar.model(m).data(d).name])); %R2013a compatible
qmodelparacond = ismember(strrep(ptmp, '_filename', ['_' ar.model(m).data(d).name]), ar.model(m).data(d).p); %R2013a compatible
ar.model(m).data(d).fp(qcondparamodel) = strrep(ar.model(m).fp(qmodelparacond), '_filename', ['_' ar.model(m).data(d).name]);

if ( substitutions == 1 )
    % Substitution code path (beta)
    from        = {};
    to          = {};
    ismodelpar  = [];
    
    % Fetch desired substitutions
    while(~isempty(C{1}) && ~strcmp(C{1},'RANDOM'))
        arFprintf( 3, '.' );
        from(end+1)         = C{1}; %#OK<AGROW>
        to(end+1)           = C{2}; %#OK<AGROW>
        ismodelpar(end+1)   = sum(ismember(ar.model(m).data(d).p, C{1})); %#OK<AGROW>
        
        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
    end    
    
    % Perform selfsubstitutions
    to = arSubsRepeated( to, fromSubs, toSubs, str2double(matVer.Version) );
    
    % Store substitutions in ar structure
    for a = 1 : length( from )
        qcondpara = ismember(ar.model(m).data(d).p, from{a}); %R2013a compatible
        if(sum(qcondpara)>0)
            ar.model(m).data(d).fp{qcondpara} = ['(' to{a} ')'];
        else
            warning('unknown parameter in conditions: %s (did you mean to place it under SUBSTITUTIONS?)', from{a}); %#ok<WNTAG>
        end
    end
else
    % old code path
    while(~isempty(C{1}) && ~strcmp(C{1},'RANDOM'))
        arFprintf( 3, '.' );
        qcondpara = ismember(ar.model(m).data(d).p, C{1}); %R2013a compatible
        if(sum(qcondpara)>0)
            ar.model(m).data(d).fp{qcondpara} = ['(' cell2mat(C{2}) ')'];
        elseif(strcmp(cell2mat(C{1}),'PARAMETERS'))
        else
            warning('unknown parameter in conditions %s', cell2mat(C{1}));
        end
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
    end
end

% extra conditional parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fp, 'UniformOutput', false);
ar.model(m).data(d).pcond = setdiff(vertcat(varlist{:}), ar.model(m).data(d).p); %R2013a compatible
      
% collect parameters conditions
pcond = union(ar.model(m).data(d).p, ar.model(m).data(d).pcond); %R2013a compatible

% RANDOM
arFprintf( 3, '[ OK ]\nReading randoms' );
if ~isempty( ar.model(m).prand )
    ar.model(m).data(d).prand = ar.model(m).prand;
    ar.model(m).data(d).rand_type = ar.model(m).rand_type;    
else
    ar.model(m).data(d).prand = {};
    ar.model(m).data(d).rand_type = [];
end
[C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
while(~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS'))
    arFprintf( 3, '.' );
    ar.model(m).data(d).prand{end+1} = cell2mat(C{1});
    if(strcmp(C{2}, 'INDEPENDENT'))
        ar.model(m).data(d).rand_type(end+1) = 0;
    elseif(strcmp(C{2}, 'NORMAL'))
        ar.model(m).data(d).rand_type(end+1) = 1;
    else
        warning('unknown random type %s', cell2mat(C{2}));  %#ok<WNTAG>
    end
    [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
end

if ( opts.expsplit )
    ar.model(m).data(d).rand_type(end+1) = 0;
    ar.model(m).data(d).prand{end+1} = opts.expsplit_args;
end

% PARAMETERS
arFprintf( 3, '[ OK ]\nReading parameters...' );
if(~isfield(ar, 'pExternLabels'))
    ar.pExternLabels = {};
    ar.pExtern = [];
    ar.qFitExtern = [];
    ar.qLog10Extern = [];
    ar.lbExtern = [];
    ar.ubExtern = [];
end
[C, fid] = arTextScan(fid, '%s %f %n %n %n %n\n',1, 'CommentStyle', ar.config.comment_string);
while(~isempty(C{1}))
    ar.pExternLabels(end+1) = C{1};
    ar.pExtern(end+1) = C{2};
    ar.qFitExtern(end+1) = C{3};
    ar.qLog10Extern(end+1) = C{4};
    ar.lbExtern(end+1) = C{5};
    ar.ubExtern(end+1) = C{6};
    [C, fid] = arTextScan(fid, '%s %f %n %n %n %n\n',1, 'CommentStyle', ar.config.comment_string);
end

% plot setup
arFprintf( 3, '[ OK ]\nPlot setup...' );
if(isfield(ar.model(m).data(d), 'response_parameter') && ...
        ~isempty(ar.model(m).data(d).response_parameter))
    if(sum(ismember(ar.model(m).data(d).p ,ar.model(m).data(d).response_parameter))==0 && ... %R2013a compatible
            sum(ismember(ar.model(m).data(d).pcond ,ar.model(m).data(d).response_parameter))==0) %R2013a compatible
        arParsingError( fid, 'invalid response parameter %s', ar.model(m).data(d).response_parameter);
    end
end
if(~isfield(ar.model(m), 'plot'))
    ar.model(m).plot(1).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
else
    ar.model(m).plot(end+1).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
end
ar.model(m).plot(end).doseresponse = ar.model(m).data(d).doseresponse;
ar.model(m).plot(end).doseresponselog10xaxis = true;
ar.model(m).plot(end).dLink = d;
ar.model(m).plot(end).ny = length(ar.model(m).data(d).y);
ar.model(m).plot(end).condition = {};
jplot = length(ar.model(m).plot);

if ( ~isstruct( fid ) )
    fclose(fid);
end

% XLS file
arFprintf( 3, 'Read def file [ OK ]\n' );
if(~strcmp(extension,'none') && ( ...
    (exist([DataPath, name '.xlsx'],'file') && strcmp(extension,'xlsx')) ||...
    (exist([DataPath, name '.xls'],'file') && strcmp(extension,'xls')) || ...
    (exist([DataPath, name '.csv'],'file') && strcmp(extension,'csv'))))
    arFprintf(2, 'loading data #%i, from file %s%s.%s...\n', d, DataPath, name, extension);
    dataFound = true;

    % read from file
    if(contains(extension,'xls'))
        warntmp = warning;
        warning('off','all')
        
        arFprintf( 3, '[ OK ]\nBegin reading data (xls) ...' );
        if (exist([DataPath, name '.xls'],'file'))      
            [data, Cstr] = xlsread([DataPath, name '.xls']);
        elseif (exist([DataPath, name '.xlsx'],'file'))      
            [data, Cstr] = xlsread([DataPath, name '.xlsx']);
        end
        arFprintf( 3, '[ OK ]\n' );
        
        if(length(data(1,:))>length(Cstr(1,:)))
            data = data(:,1:length(Cstr(1,:)));
        end
        
        warning(warntmp);
        
        timevar = Cstr(1,1);
        header = Cstr(1,2:end);
        header = strrep(header,' ',''); % remove spaces which are sometimes in the column header by accident    
        times = data(:,1);
        qtimesnonnan = ~isnan(times);
        times = times(qtimesnonnan);
        data = data(qtimesnonnan,2:end);
        if(size(data,2)<length(header))
            data = [data nan(size(data,1),length(header)-size(data,2))];
        end
        
        Cstr = Cstr(2:end,2:end);
        dataCell = cell(size(data));
        for j1 = 1:size(data,1)
            for j2 = 1:size(data,2)
                if(isnan(data(j1,j2)))
                    if(j1<=size(Cstr,1) && j2<=size(Cstr,2) && ~isempty(Cstr{j1,j2}))
                        dataCell{j1,j2} = Cstr{j1,j2};
                    else
                        dataCell{j1,j2} = header{j2};
                    end
                else
                    dataCell{j1,j2} = num2str(data(j1,j2));
                end
            end
        end
        
    elseif(strcmp(extension,'csv'))
        arFprintf( 3, '[ OK ]\nBegin reading data (csv) ...' );
        [header, data, dataCell] = arReadCSVHeaderFile([DataPath, name '.csv'], ',', true);
        arFprintf( 3, '[ OK ]\n' );
        
        timevar = strtrim(header(1));
        header = header(2:end);
        times = data(:,1);
        data = data(:,2:end);
        dataCell = dataCell(:,2:end);
    end
    
    % remove time points that we don't want
    if ( opts.removeconditions )
        arFprintf( 3, 'Removing undesired time points...' );
        selected = true(1, size(times,1));
        if ( opts.removeconditions )
            for a = 1 : 2 : length( opts.removeconditions_args )
                if ( strcmp( timevar, opts.removeconditions_args{a} ) )
                    % If the argument is a function handle, we evaluate them
                    % for each element
                    val = opts.removeconditions_args{a+1};
                    if ( isa(val, 'function_handle') )
                        for jv = 1 : length( times )
                            accepted(jv) = val(num2str(times(jv)));
                        end
                    else
                        arParsingError( fid,  'Filter argument for removecondition is of the wrong type' );
                    end
                    selected = selected & ~accepted;
                end
            end
        end
        times = times(selected);
        data = data(selected,:);
        dataCell = dataCell(selected,:);
    end
      
    % random effects
    arFprintf( 3, 'Processing random effects...' );
    prand = ar.model(m).data(d).prand;
    if(opts.splitconditions)
        prand = union(prand, opts.splitconditions_args);
    end
    qrandis = ismember(header, prand); %R2013a compatible
    if(sum(qrandis) > 0)
        qobs = ismember(header, ar.model(m).data(d).y); %R2013a compatible
        
        randis_header = header(qrandis);
        qrandis_header_nosplit = ismember(randis_header, ar.model(m).data(d).prand);
        
        if ~isempty(dataCell)
            [randis, ~, jrandis] = uniqueRowsCA(dataCell(:,qrandis));
        else
            [randis, ~, jrandis] = unique(data(:,qrandis),'rows');
            randis = cellstr(num2str(randis));
        end
               
        for j=1:size(randis,1)
            qvals = jrandis == j;
            tmpdata = data(qvals,qobs);
            if(sum(~isnan(tmpdata(:)))>0 || ~removeEmptyObs)
                arFprintf(2, 'local random effect #%i:\n', j)
                
                if(j < size(randis,1))
                    ar.model(m).data(d+1) = ar.model(m).data(d);
                    ar.model(m).plot(jplot+1) = ar.model(m).plot(jplot);
                end
                
                pcondmod = pcond;
                for jj=1:size(randis,2)
                    if(qrandis_header_nosplit(jj))
                        arFprintf(2, '\t%20s = %s\n', randis_header{jj}, randis{j,jj})
                        
                        ar.model(m).plot(jplot).name = [ar.model(m).plot(jplot).name '_' ...
                            randis_header{jj} randis{j,jj}];
                        
                        ar.model(m).data(d).name = [ar.model(m).data(d).name '_' ...
                            randis_header{jj} randis{j,jj}];
                        ar.model(m).data(d).fprand = randis{j,jj};
                        
                        ar.model(m).data(d).fy = strrep(ar.model(m).data(d).fy, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).py = strrep(ar.model(m).data(d).py, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        
                        ar.model(m).data(d).fystd = strrep(ar.model(m).data(d).fystd, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).pystd = strrep(ar.model(m).data(d).pystd, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        
                        %                     ar.model(m).data(d).p = strrep(ar.model(m).data(d).p, ...
                        %                         randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).fp = strrep(ar.model(m).data(d).fp, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).pcond = strrep(ar.model(m).data(d).pcond, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        
                        for jjj=1:length(ar.model(m).data(d).py_sep)
                            ar.model(m).data(d).py_sep(jjj).pars = strrep(ar.model(m).data(d).py_sep(jjj).pars, ...
                                randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        end
                        
                        pcondmod = strrep(pcondmod, randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    else
                        arFprintf(2, '\t%20s (split only)\n', randis_header{jj})
                        
                        ar.model(m).plot(jplot).name = [ar.model(m).plot(jplot).name '_' ...
                            randis_header{jj} randis{j,jj}];
                    end
                end
                
                arFprintf( 4, 'Setting conditions ...\n' );
                if ~isempty(dataCell)
                    [ar,d,fail] = setConditions(fid, ar, m, d, jplot, header, times(qvals), data(qvals,:), dataCell(qvals,:), ...
                        pcondmod, removeEmptyObs, dpPerShoot, opts);
                else
                    [ar,d,fail] = setConditions(fid, ar, m, d, jplot, header, times(qvals), data(qvals,:), dataCell, ...
                        pcondmod, removeEmptyObs, dpPerShoot, opts);
                end
                arFprintf( 4, 'Condition set ... [ OK ]\n' );
                
                % Only increment if some data was actually set.
                if (~fail)
                    if(j < size(randis,1))
                        d = d + 1;
                        jplot = jplot + 1;
                        ar.model(m).plot(jplot).dLink = d;
                    end
                
                    % Check whether the user specified any variables with reserved words.
                    checkReserved(m, d);
                else
                    % Remove file which failed to provide any data
                    fprintf(2, 'local random effect #%i: no matching data (%d), removed\n', j ,d);
                    ar.model.data(d) = [];
                end
                
            else
                arFprintf(2, 'local random effect #%i: no matching data, skipped\n', j);
            end
        end
    else
        arFprintf( 4, 'Setting conditions ...\n' );
        ar = setConditions(fid, ar, m, d, jplot, header, times, data, dataCell, pcond, removeEmptyObs, dpPerShoot, opts);
        arFprintf( 4, 'Condition set ... [ OK ]\n' );
        
        % Check whether the user specified any variables with reserved words.
        checkReserved(m, d);
    end
else
    dataFound = false;
    warning('Cannot find data file corresponding to %s', name);
    ar.model(m).data(d).condition = [];
end


% remember the function call
ar.setup.commands{end+1} = mfilename; % this file name
ar.setup.arguments{end+1} = {name,m,extension, removeEmptyObs, varargin{:}}; % 
if dataFound
    ar.setup.datafiles{end+1} = {[DataPath,name,'.def'],[DataPath,name,'.',extension]};
else
    ar.setup.datafiles{end+1} = {[DataPath,name,'.def'],''};
end
ar.setup.modelfiles{end+1} = '';

% sort fields
ar = orderfields(ar);
ar.model = orderfields(ar.model);
ar.model(m).data = orderfields(ar.model(m).data);
ar.model(m).plot = orderfields(ar.model(m).plot);





function checkReserved(m, d)
    global ar;

    % Check whether the user specified any variables with reserved words.
    for a = 1 : length( ar.model(m).data(d).fu )
        arCheckReservedWords( symvar(ar.model(m).data(d).fu{a}), sprintf( 'input function of %s', ar.model(m).data(d).name ), ar.model(m).u{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fy )
        arCheckReservedWords( symvar(ar.model(m).data(d).fy{a}), sprintf( 'observation function of %s', ar.model(m).data(d).name ), ar.model(m).data(d).y{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fystd )
        arCheckReservedWords( symvar(ar.model(m).data(d).fystd{a}), sprintf( 'observation standard deviation function of %s', ar.model(m).data(d).name ), ar.model(m).data(d).y{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fp )
        arCheckReservedWords( symvar(ar.model(m).data(d).fp{a}), sprintf( 'condition parameter transformations of %s', ar.model(m).data(d).name ), ar.model(m).data(d).p{a} );
    end   
    arCheckReservedWords( ar.model(m).data(d).p, 'parameters' );
    arCheckReservedWords( ar.model(m).data(d).y, 'observable names' );

function [ar,d, fail] = setConditions(fid, ar, m, d, jplot, header, times, data, dataCell, pcond, removeEmptyObs, dpPerShoot, opts)

% matVer = ver('MATLAB');

% normalization of columns
fail = 0;
nfactor = max(data, [], 1);

qobs = ismember(header, ar.model(m).data(d).y) & sum(~isnan(data),1)>0; %R2013a compatible
qhasdata = ismember(ar.model(m).data(d).y, header(qobs)); %R2013a compatible

% conditions
if (~opts.removeconditions)
    qcond = ismember(header, pcond); %R2013a compatible
else
    % Add the condi's we force filtering over (override)
    qcond = ismember(header, pcond) | ismember(header, opts.removeconditions_args(1:2:end)); %R2013a compatible
end

% Refine dose responses if requested
if ( opts.resampledoseresponse )
    resolution = opts.resamplingresolution;
    if ( ar.model(m).data(d).doseresponse == true )
        responsePar = ismember( header, ar.model(m).data(d).response_parameter );
        if ( sum( responsePar ) )
            fprintf( '  => Refining dose response for %s\n', ar.model(m).data(d).response_parameter );
            if ( sum( responsePar ) > 1 )
                arParsingError( fid,  'Response parameter ambiguous during dose response refinement' );
            end

            % Which columns define the conditions
            conds = qcond & ~responsePar;

            % Grab unique conditions (note the inclusion of time)
            [uniqueCondi, ~, ib] = unique( [times data( :, qcond & ~responsePar ) ], 'rows' );
            nConditions = size(uniqueCondi, 1);

            % For each unique condition determine the maximum and minimum value
            % of the response parameter
            extraData = NaN(nConditions*resolution, size(data,2));
            extraTimes = NaN(nConditions*resolution, 1);

            for jui = 1 : size( uniqueCondi, 1 )
                dataChunk = data( ib == jui, responsePar );
                mi = min( dataChunk );
                ma = max( dataChunk );
                extraPoints = mi : (ma-mi)/(resolution-1) : ma;

                if ( opts.refinelog )     
                    mi = log10( max( [mi, 1e-8] ) );
                    ma = log10( max( [ma, 1e-8] ) );
                    extraPoints = 10.^(mi : (ma-mi)/(resolution-1) : ma);
                end
                
                % Fix to make sure the length is correct for the next
                % assignment. This is ok, since the unique will remove
                % these points again at a later stage.
                if ( mi == ma )
                    extraPoints = repmat( mi, 1, resolution );
                end

                % Fill extra data with current condition
                condition = uniqueCondi(jui,:);
                extraData((jui-1)*resolution+1:jui*resolution, conds) = repmat(condition(2:end), resolution, 1);
                extraTimes((jui-1)*resolution+1:jui*resolution) = repmat(condition(1), resolution, 1);

                % Fill the dependent variable with the new dependent variable values
                extraData((jui-1)*resolution+1:jui*resolution, responsePar) = extraPoints;
            end
        end
        data = [ data ; extraData ];
        dataCell = [ dataCell; cellfun(@num2str,num2cell(extraData), 'UniformOutput', false) ];
        times = [ times ; extraTimes ];
    end
end

if(sum(qcond) > 0)
    condi_header = header(qcond);
    if ~isempty(dataCell)
        [condis, ind, jcondis] = uniqueRowsCA(dataCell(:,qcond));
    else
        [condis, ind, jcondis] = unique(data(:,qcond),'rows');
        condis = mymat2cell(condis);
    end

    if (opts.removeconditions || opts.removeemptyconds)
        selected = true(1, size(condis,1));
        if ( opts.removeconditions )
            for a = 1 : 2 : length( opts.removeconditions_args )
                cc = ismember( condi_header, opts.removeconditions_args{a} );
                if ( sum( cc ) > 0 )
                    values = condis(:,cc);

                    % If the argument is a function handle, we evaluate them
                    % for each element
                    val = opts.removeconditions_args{a+1};
                    if ( isa(val, 'function_handle') )
                        for jv = 1 : length( values )
                            accepted(jv) = val(values{jv});
                        end
                    else
                        if (isnumeric(val))
                            val = num2str(val);
                        end
                        if ~ischar(val)
                            arParsingError( fid,  'Filter argument for removecondition is of the wrong type' );
                        end
                        accepted = ismember(values, val).';
                    end
                    selected = selected & ~accepted;
                end
            end
        end
        if(opts.removeemptyconds)
            % Find out for which conditions we actually have data
            hasD = max(~isnan(data(ind,qobs)), [], 2);
            selected(hasD==0) = false;
        end
        condis = condis(selected,:);
        
        % Recompute jcondi's (list which points which data row corresponds
        % to which condition.
        mapTo   = cumsum(selected);
        mapTo(~selected) = -1;
        jcondis = mapTo(jcondis);
    end
    
    % exit if no data left
    if(size(condis,1)==0)
        fail = 1;
        return
    end
       
    active_condi = false(size(condis(1,:)));
    tmpcondi = condis(1,:);
    for j1=2:size(condis,1)
        for j2=1:size(condis,2)
            active_condi(j2) = active_condi(j2) | (~strcmp(tmpcondi{j2}, condis{j1,j2}));
        end
    end
        
    for j=1:size(condis,1)
        
        arFprintf(2, 'local condition #%i:\n', j)
        
        if(j < size(condis,1))
            if(length(ar.model(m).data) > d)
                ar.model(m).data(d+2) = ar.model(m).data(d+1);
            end
            ar.model(m).data(d+1) = ar.model(m).data(d);
        end
        
        % remove obs without data
        if(removeEmptyObs)
            for jj=find(~qhasdata)
                arFprintf(2, '\t%20s no data, removed\n', ar.model(m).data(d).y{jj});
                jjjs = find(ismember(ar.model(m).data(d).p, ar.model(m).data(d).py_sep(jj).pars)); %R2013a compatible
                jjjs = jjjs(:)';
                for jjj=jjjs
                    remove = 1;
                    for jjjj = find(qhasdata)
                        if sum(ismember(ar.model(m).data(d).py_sep(jjjj).pars, ar.model(m).data(d).p(jjj))) > 0 %R2013a compatible
                            remove = 0;
                        end
                    end
                    if remove
                        ar.model(m).data(d).fp{jjj} = '0';
                    end
                end
            end
            ar.model(m).data(d).y = ar.model(m).data(d).y(qhasdata);
            ar.model(m).data(d).yNames = ar.model(m).data(d).yNames(qhasdata);
            ar.model(m).data(d).yUnits = ar.model(m).data(d).yUnits(qhasdata,:);
            ar.model(m).data(d).normalize = ar.model(m).data(d).normalize(qhasdata);
            ar.model(m).data(d).logfitting = ar.model(m).data(d).logfitting(qhasdata);
            ar.model(m).data(d).logplotting = ar.model(m).data(d).logplotting(qhasdata);
            ar.model(m).data(d).fy = ar.model(m).data(d).fy(qhasdata);
            ar.model(m).data(d).fystd = ar.model(m).data(d).fystd(qhasdata);
            ar.model(m).data(d).py_sep = ar.model(m).data(d).py_sep(qhasdata);
        end
        
        for jj=1:size(condis,2)
            if(~isempty(condis{j,jj}))
                arFprintf(2, '\t%20s = %s\n', condi_header{jj}, condis{j,jj})
                
                qcondjj = ismember(ar.model(m).data(d).p, condi_header{jj}); %R2013a compatible
                if(sum(qcondjj)>0)
                    ar.model(m).data(d).fp{qcondjj} =  ['(' condis{j,jj} ')'];
                end
                qcondjj = ~strcmp(ar.model(m).data(d).p, ar.model(m).data(d).fp');
                if(~isnan(str2double(condis{j,jj})))
%                     ar.model(m).data(d).fp(qcondjj) = strrep(ar.model(m).data(d).fp(qcondjj), ...
%                         condi_header{jj}, condis{j,jj});

                    ar.model(m).data(d).fp(qcondjj) = regexprep(ar.model(m).data(d).fp(qcondjj),...
                        sprintf('\\<%s\\>', condi_header{jj}),condis{j,jj});
                    
%                     tmpfp = subs(sym(ar.model(m).data(d).fp(qcondjj)), ...
%                         sym(condi_header{jj}), sym(condis{j,jj}));
%                     jps = find(qcondjj);
%                     for jp = 1:length(jps)
%                         ar.model(m).data(d).fp{jps(jp)} = char(tmpfp(jp));
%                     end
                end
                
                ar.model(m).data(d).condition(jj).parameter = condi_header{jj};
                ar.model(m).data(d).condition(jj).value = condis{j,jj};
                
                % plot
                if(active_condi(jj))
                    if(ar.model(m).data(d).doseresponse==0 || ~strcmp(condi_header{jj}, ar.model(m).data(d).response_parameter))
                        if(length(ar.model(m).plot(jplot).condition) >= j && ~isempty(ar.model(m).plot(jplot).condition{j}))
                            ar.model(m).plot(jplot).condition{j} = [ar.model(m).plot(jplot).condition{j} ' & ' ...
                                ar.model(m).data(d).condition(jj).parameter '=' ...
                                ar.model(m).data(d).condition(jj).value];
                        else
                            ar.model(m).plot(jplot).condition{j} = [ar.model(m).data(d).condition(jj).parameter '=' ...
                                ar.model(m).data(d).condition(jj).value];
                        end
                    end
                end
            end
        end
        
        qvals = jcondis == j;
        ar = setValues(fid, ar, m, d, header, nfactor, data(qvals,:), times(qvals));
        ar.model(m).data(d).tLim(2) = round(max(times)*1.1);
        
        if(dpPerShoot~=0)
            [ar,d] = doMS(ar,m,d,jplot,dpPerShoot);
        end
        
        if(j < size(condis,1))
            d = d + 1;
            ar.model(m).plot(jplot).dLink(end+1) = d;
        end
    end
else
    ar.model(m).data(d).condition = [];
    
    % remove obs without data
    if(removeEmptyObs)
        for jj=find(~qhasdata)
            arFprintf(2, '\t%20s no data, removed\n', ar.model(m).data(d).y{jj});
            jjjs = find(ismember(ar.model(m).data(d).p, ar.model(m).data(d).py_sep(jj).pars)); %R2013a compatible
            jjjs = jjjs(:)';
            for jjj=jjjs
                remove = 1;
                for jjjj = find(qhasdata)
                    if sum(ismember(ar.model(m).data(d).py_sep(jjjj).pars, ar.model(m).data(d).p(jjj))) > 0 %R2013a compatible
                        remove = 0;
                    end
                end
                if(remove==1)
                    ar.model(m).data(d).fp{jjj} = '0';
                end
            end
        end
        ar.model(m).data(d).y = ar.model(m).data(d).y(qhasdata);
        ar.model(m).data(d).yNames = ar.model(m).data(d).yNames(qhasdata);
        ar.model(m).data(d).yUnits = ar.model(m).data(d).yUnits(qhasdata,:);
        ar.model(m).data(d).normalize = ar.model(m).data(d).normalize(qhasdata);
        ar.model(m).data(d).logfitting = ar.model(m).data(d).logfitting(qhasdata);
        ar.model(m).data(d).logplotting = ar.model(m).data(d).logplotting(qhasdata);
        ar.model(m).data(d).fy = ar.model(m).data(d).fy(qhasdata);
        ar.model(m).data(d).fystd = ar.model(m).data(d).fystd(qhasdata);
        ar.model(m).data(d).py_sep = ar.model(m).data(d).py_sep(qhasdata);
    end
    
    ar = setValues(fid, ar, m, d, header, nfactor, data, times);
    ar.model(m).data(d).tLim(2) = round(max(times)*1.1);
    
    if(dpPerShoot~=0)
        [ar,d] = doMS(ar,m,d,jplot,dpPerShoot);
    end
end




function C = mymat2cell(D)
C = cell(size(D));
for j=1:size(D,1)
    for jj=1:size(D,2)
        C{j,jj} = num2str(D(j,jj));
    end
end

function [ar,d] = doMS(ar,m,d,jplot,dpPerShoot)

tExp = ar.model(m).data(d).tExp;

if(dpPerShoot ~= 1)
    nints = ceil(length(tExp) / dpPerShoot);
    tboarders = linspace(min(tExp),max(tExp),nints+1);
else
    tboarders = union(0,tExp); %R2013a compatible
    nints = length(tboarders)-1;
end

if(nints==1)
    return;
end

arFprintf(2, 'using %i shooting intervals\n', nints);
ar.model(m).ms_count = ar.model(m).ms_count + 1;
ar.model(m).data(d).ms_index = ar.model(m).ms_count;

for j=1:nints
    ar.model(m).data(d).ms_snip_index = j;
    if(j<nints)
        ar.model(m).data(end+1) = ar.model(m).data(d);
        ar.model(m).plot(jplot).dLink(end+1) = d+1;
    end
    
    if(j>1)
        ar.ms_count_snips = ar.ms_count_snips + 1;       
        qtodo = ismember(ar.model(m).data(d).p, ar.model(m).px0); %R2013a compatible
        ar.model(m).data(d).fp(qtodo) = strrep(ar.model(m).data(d).p(qtodo), 'init_', sprintf('init_MS%i_', ar.ms_count_snips));
    end
    
    if(j<nints)
        ar.model(m).data(d).tExp = ar.model(m).data(d).tExp(tExp>=tboarders(j) & tExp<tboarders(j+1));
        ar.model(m).data(d).yExp = ar.model(m).data(d).yExp(tExp>=tboarders(j) & tExp<tboarders(j+1),:);
        ar.model(m).data(d).yExpStd = ar.model(m).data(d).yExpStd(tExp>=tboarders(j) & tExp<tboarders(j+1),:);
    else
        ar.model(m).data(d).tExp = ar.model(m).data(d).tExp(tExp>=tboarders(j) & tExp<=tboarders(j+1));
        ar.model(m).data(d).yExp = ar.model(m).data(d).yExp(tExp>=tboarders(j) & tExp<=tboarders(j+1),:);
        ar.model(m).data(d).yExpStd = ar.model(m).data(d).yExpStd(tExp>=tboarders(j) & tExp<=tboarders(j+1),:);
    end
    
    ar.model(m).data(d).tLim = [tboarders(j) tboarders(j+1)];
    ar.model(m).data(d).tLimExp = ar.model(m).data(d).tLim;
    
    if(j<nints)
        d = d + 1;
    end
end


function ar = setValues(fid, ar, m, d, header, nfactor, data, times)
ar.model(m).data(d).tExp = times;
ar.model(m).data(d).yExp = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStd = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpRaw = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStdRaw = nan(length(times), length(ar.model(m).data(d).y));

for j=1:length(ar.model(m).data(d).y)
    q = ismember(header, ar.model(m).data(d).y{j}); %R2013a compatible
    
    if(sum(q)==1)
        ar.model(m).data(d).yExp(:,j) = data(:,q);
        ar.model(m).data(d).yExpRaw(:,j) = data(:,q);
        arFprintf(2, '\t%20s -> %4i data-points assigned', ar.model(m).data(d).y{j}, sum(~isnan(data(:,q))));
        
        % normalize data
        if(ar.model(m).data(d).normalize(j))
            ar.model(m).data(d).yExp(:,j) = ar.model(m).data(d).yExp(:,j) / nfactor(q);
            arFprintf(2, ' normalized');
        end
        
        % log-fitting
        if(ar.model(m).data(d).logfitting(j))
            qdatapos = ar.model(m).data(d).yExp(:,j)>0;
            nancount = length(isnan(ar.model(m).data(d).yExp(:,j)));
            ar.model(m).data(d).yExp(qdatapos,j) = log10(ar.model(m).data(d).yExp(qdatapos,j));
            ar.model(m).data(d).yExp(~qdatapos,j) = nan;
            if(sum(~qdatapos)==0)
                arFprintf(2, ' for log-fitting');
            else
                arFprintf(2, ' for log-fitting (%i values <=0 removed, %i NaN values removed)', sum(~qdatapos)-nancount, nancount);
            end
        end
        
        % empirical stds
        qstd = ismember(header, [ar.model(m).data(d).y{j} '_std']); %R2013a compatible
        if(sum(qstd)==1)
            ar.model(m).data(d).yExpStdRaw(:,j) = data(:,qstd);
            ar.model(m).data(d).yExpStd(:,j) = data(:,qstd);
            arFprintf(2, ' with stds');
            if(ar.model(m).data(d).normalize(j))
                ar.model(m).data(d).yExpStd(:,j) = ar.model(m).data(d).yExpStd(:,j) / nfactor(q);
                arFprintf(2, ' normalized');
            end
        elseif(sum(qstd)>1)
            arParsingError( fid, 'multiple std colums for observable %s', ar.model(m).data(d).y{j})
        end
        
    elseif(sum(q)==0)
        arFprintf(2, '*\t%20s -> not assigned', ar.model(m).data(d).y{j});
    else
        arParsingError( fid, 'multiple data colums for observable %s', ar.model(m).data(d).y{j})
    end
    
    arFprintf(1, '\n');
end

function num = checkNum( num, defaultValue )
    if ( ~isnumeric( num ) || isempty( num ) || isnan( num ) )
        num = defaultValue;
    end