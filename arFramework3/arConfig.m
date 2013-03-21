function varargout = arConfig(varargin)
% ARTUNER2 M-file for arConfig.fig
%      ARCONFIG, by itself, creates a new ARCONFIG or raises the existing
%      singleton*.
%
%      H = ARCONFIG returns the handle to a new ARCONFIG or the handle to
%      the existing singleton*.
%
%      ARCONFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARCONFIG.M with the given input arguments.
%
%      ARCONFIG('Property','Value',...) creates a new ARCONFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arConfig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arConfig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arConfig

% Last Modified by GUIDE v2.5 17-Aug-2010 17:10:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arConfig_OpeningFcn, ...
                   'gui_OutputFcn',  @arConfig_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before arConfig is made visible.
function arConfig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to arConfig (see VARARGIN)

% Choose default command line output for arConfig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes arConfig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = arConfig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_table(hObject)
global ar
data = cell(length(ar.p), 6);
for j=1:length(ar.p)
    data{j,1} = ar.pLabel{j};
    data{j,2} = ar.p(j);
    if(ar.qFit(j)==0)
        data{j,3} = 'fixed';
    elseif(ar.qFit(j)==1)
        data{j,3} = 'free';
    elseif(ar.qFit(j)==2)
        data{j,3} = 'constant';
    end
    data{j,4} = ar.qLog10(j)==1;
    data{j,5} = ar.lb(j);
    data{j,6} = ar.ub(j);
end
set(hObject, 'Data', data);


% --- Executes during object creation, after setting all properties.
function uitable2_CreateFcn(hObject, eventdata, handles)
colnames = {'parameter', 'value', 'fit', 'log10', 'lb', 'ub'};
colwidth = {200 80 85 35 35 35};
columnformat = {'char', 'numeric', {'free' 'fixed', 'constant'}, 'logical', 'numeric', 'numeric'}; 
columneditable =  [false true true true true true]; 

set(hObject, 'ColumnName', colnames);
set(hObject, 'ColumnWidth', colwidth);
set(hObject, 'ColumnEditable', columneditable)
set(hObject, 'ColumnFormat', columnformat);

update_table(hObject);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
arFit(true);
update_table(handles.uitable2);
arChi2(false)
arPlot(false, true)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
arFitObs(true);
update_table(handles.uitable2);
arChi2(false)
arPlot(false, true)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
arFitDyn(true);
update_table(handles.uitable2);
arChi2(false)
arPlot(false, true)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global ar
if(isfield(ar, 'tuner') && isfield(ar.tuner, 'index2') && ~isempty(ar.tuner.index2))
    arFitSingle(ar.tuner.index2, true);
    update_table(handles.uitable2);
    arChi2(false)
    arPlot(false, true)
end

% --- Executes when selected cell(s) is changed in uitable2.
function uitable2_CellSelectionCallback(hObject, eventdata, handles)
global ar
if(~isempty(eventdata.Indices) && eventdata.Indices(2)==1)
    ar.tuner.index2 = eventdata.Indices(1);
else
    ar.tuner.index2 = [];
end

% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global ar
if(~isempty(eventdata.Indices) && eventdata.Indices(2)==2)
    if(~isnan(eventdata.NewData))
        if(eventdata.NewData < ar.lb(ar.tuner.index))
            eventdata.NewData = ar.lb(ar.tuner.index);
        end
        if(eventdata.NewData > ar.ub(ar.tuner.index))
            eventdata.NewData = ar.ub(ar.tuner.index);
        end
        ar.p(eventdata.Indices(1)) = eventdata.NewData;
        % arChi2(false)
        arPlot(false, true)
    end
end
