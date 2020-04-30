
% arLoadData(name, [m], [extension], [removeEmptyObs], [mixingCond], [fieldsvar], [opts])
%
% Load data set to next free slot(s)
%
%   name                      filename of data definition file
%   m                         target position (int) for modelor: model 
%                             name (string)
%                             [last loaded model]
%   extension                 data file name-extension: 'xls', 'csv',         
%                             'none' = don't load data
%                             ['xls']
%   mixingCond       double   mixing condition = [x x x]
%                             determines how the user defined data struct
%                             will be added to the dada def read from file
%                             elemenst are either zero or on
%                             1st element: data struct
%                                   0: discars new data stuct
%                                      other elements are not important
%                                   1: takes into account the data struct  
%                             2nd element: def information
%                                   0: updates the data struct by the
%                                      new one
%                                   1: replace the existing data struct
%                                          with the new one 
%                             3rd element: data file
%                                   0: updates the data file by the
%                                      new one
%                                   1: replace the existing data file
%                                      with the new one
%                             [0 0 0]
%                             Short form: 0=[0 0 0], 1=[1 1 1]
%   removeEmptyObs            remove observation without data
%                             [false]
%   fieldsvar                 data structure created by user. it is defined exactly 
%                             like you are building def file with the same sections.
%                             for more information check arCreateDataStructure.       
%
%   opts                      additional option flags
%
%   Optional option flags are:
%   'RemoveConditions'        This flag followed by a list of conditions will
%                             allow you to filter the data that you load. Note
%                             that the function takes both values, strings
%                             or function handles. When you provide a function
%                             handle, the function will be evaluated for each
%                             condition value (in this case input_dcf). You
%                             should make the function return 1 if the condition
%                             is to be removed. Note that the input to the
%                             function is a *string* not a number.
%                             Example:
%                               arLoadData( 'mydata', 1, 'csv', true, 'RemoveConditions', ...
%                               {'input_il6', '0', 'input_dcf', @(dcf)str2num(dcf)>0};
%   'RemoveEmptyConds'        This flag allows you to remove conditions that have
%                             no data points.
%   'ResampleDoseResponse'    This flag allows you to increase the resolution
%                             of dose responses for plotting purposes. Note
%                             that this makes model evaluation and
%                             compilation much more time intensive (hence it
%                             is off by default).
%   'ResamplingResolution'    Number of extra points to interpolate dose
%                             response (default: 25)
%   'RefineLog'               Perform the refinement on a log scale
%   'expsplit'                Split replicate specific parameters based on the 
%                             column identifier specified in the next argument
%                             Example:
%                               arLoadData('mRNA/mRNA_pretreatment', 1, 'csv', true, 'expsplit', 'nExpID');
%                               will recognize parameters with nExpID in the name 
%                               and replace them with nExpID0, nExpID1 etc. if
%                               the excel data file has a column named nExpID  
%                               which specifies these values.
%
%   'RemoveObservables'       This can be used to omit observables. Simply pass a
%                             cell array with names of observables that should be 
%                             ignored in the data.
%
%   'IgnoreInputs'            Ignore input overrides. Simply pass a cell array
%                             with names of inputs that should be ignored (the 
%                             model default will be used instead).
% 
%   'DataPath'                Path to the data files.
%                             Default: DataPath = 'Data/'
% 
%   The data file specification is as follows:
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
%   Defining the measurement noise in the data sheet:
%   Instead of using a parametrized function for the measurement noise, 
%   the amount of noise can be put value in the data sheet with column 
%   header same as the observable name but with the addition _std 
%   (e.g. tSTAT5_au_std and pSTAT5_au_std). These values will be 
%   interpreted as standard deviation of the data not as variance!
%
%   The following flags controll how the experimental noise enters 
%   calculations and how it is plotted:
%   ar.config.fiterror     0: Data is plotted as fitted (ar.config.fiterors).
%                         -1: Plot prediction bands as calculated by PPL.
%                         -2: Do not plot errors
%   ar.config.fiterrors   -1: Only experimental error bars used for fitting.
%                          0: Use experimental errors by default and revert
%                             to the error model for data points that have no
%                             experimental error specified.
%                          1: Only the error model is used for fitting.
%                             Experimental errors specified in the data sheet
%                             are ignored.
%   ar.config.ploterrors   0: Observables are plotted as fitted.
%                          1: Data uncertainty is plotted as error bar.
%                          2: Only error bands are plotted.   
%
% See wiki page about loading data: https://github.com/Data2Dynamics/d2d/wiki/Setting%20up%20models#2-data-definition-file
%
% See also arInit, arLoadModel 


function arLoadData(name, m, extension, removeEmptyObs, mixingCond, fieldsvar, varargin)


global ar


% Initiate parameters
if ~exist('m','var')
    m = [];
end

if ~exist('extension','var') || ~ischar(extension)
    extension='';
end

if ~exist('varargin','var') || isempty(varargin)
    opts = option();
else
    opts = option(varargin);
end


if(~exist('removeEmptyObs','var'))
    removeEmptyObs = false;
end
    
if ~exist('fieldsvar','var') || isempty(fieldsvar)
    fieldsvar={};
end


if ~exist('mixingCond','var') || isempty(mixingCond)
    mixingCond = zeros(1,3);
else
    mixingCond = setCond(mixingCond);
end


% Loads data file
RawData = arReadDataFile(name, extension, opts);


% Loads data struct from def file if requested
if ~all(mixingCond(1:2))
    DefFile = arReadDefFile(name, opts, ar.config.comment_string);
else
    DefFile={};
end
    

% Creates data structure defined explicitly by the user (without using data def file)
if mixingCond(1)==1
    CreateData = arCreateDataStructure(fieldsvar, mixingCond, name, opts);
else
    CreateData={};
end


% combines data structures and links them to data
DataDefInfo = arAddDataDefInfo(RawData, DefFile, CreateData, mixingCond, opts);


% Computes data struct and add to the global ar
if ~exist('varargin','var') || isempty(varargin)
    arAddDataStructure(name, DataDefInfo, m, removeEmptyObs, opts);
else
    arAddDataStructure(name, DataDefInfo, m, removeEmptyObs, opts, varargin);
end


% remember the function call
ar.setup.commands{end+1} = mfilename; % this file name

end


function mixingCond = setCond(mixingCond)

if  ~isnumeric(mixingCond)
    error('mixingCond should be a numeric array with three elements')
end

if length(mixingCond)>3
    error('length of input argument ''mixingCond'' must be three.')
elseif any(mixingCond>1) || any(mixingCond<0)
    error('elements of ''mixingCond'' must be either one or zero.')
    
elseif length(mixingCond)~=3
    if mixingCond(1)==0
        mixingCond(2:3)=0;
    elseif length(mixingCond)==1 && mixingCond(1)==1
        mixingCond=ones(1,3);
    else
        error('length of input argument ''mixingCond'' must be three.')
    end
end

end


function opts = option(varargin) 

switches = { 'dppershoot',       'removeconditions', 'removeobservables',    'splitconditions',...
             'removeemptyconds', 'expsplit',         'resampledoseresponse', 'resamplingresolution',...
             'refinelog',        'ignoreinputs',     'detectionlimit',       'datapath'};

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

end
