function varargout = arParameters(varargin)
% ARPARAMETERS MATLAB code for arParameters.fig
%      ARPARAMETERS, by itself, creates a new ARPARAMETERS or raises the existing
%      singleton*.
%
%      H = ARPARAMETERS returns the handle to a new ARPARAMETERS or the handle to
%      the existing singleton*.
%
%      ARPARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARPARAMETERS.M with the given input arguments.
%
%      ARPARAMETERS('Property','Value',...) creates a new ARPARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arParameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arParameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arParameters

% Last Modified by GUIDE v2.5 05-Apr-2015 19:38:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arParameters_OpeningFcn, ...
                   'gui_OutputFcn',  @arParameters_OutputFcn, ...
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


% --- Executes just before arParameters is made visible.
function arParameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to arParameters (see VARARGIN)

% Choose default command line output for arParameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes arParameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = arParameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global ar

type = {'uniform','normal','uniform with normal bounds','L1'};
C = cell(length(ar.p),9);
for jp=1:length(ar.p)
    C{jp,1} = ar.pLabel{jp};
    C{jp,2} = ar.lb(jp);
    C{jp,3} = ar.p(jp);
    C{jp,4} = ar.ub(jp);
    C{jp,5} = ar.qLog10(jp)==1;
    C{jp,6} = ar.qFit(jp)==1;
    C{jp,7} = type{ar.type(jp)+1};
    C{jp,8} = ar.mean(jp);
    C{jp,9} = ar.std(jp);
end

set(hObject, 'Data', C);


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

global ar

jp = eventdata.Indices(1);
jt = eventdata.Indices(2);
newval = eventdata.EditData;
C = get(hObject, 'Data');

type = {'uniform','normal','uniform with normal bounds','L1'};

% check for nan
if(jt~=5 && jt~=6 && jt~=7 && isnan(str2double(newval)))
    C{jp,jt} = eventdata.PreviousData;
    set(hObject, 'Data', C);
    return;
end

if(jt==2) % lb
    newval = str2double(newval);
    ar.lb(jp) = newval;
    if(ar.p(jp) < ar.lb(jp))
        ar.p(jp) = ar.lb(jp);
        C{jp,3} = ar.p(jp);
        set(hObject, 'Data', C);
    end
elseif(jt==3) % p
    newval = str2double(newval);
    if(newval < ar.lb(jp))
        newval = ar.lb(jp);
        C{jp,3} = newval;
        set(hObject, 'Data', C);
    end
    if(newval > ar.ub(jp))
        newval = ar.ub(jp);
        C{jp,3} = newval;
        set(hObject, 'Data', C);
    end
    ar.p(jp) = newval;
elseif(jt==4) % ub
    newval = str2double(newval);
    ar.ub(jp) = newval;
    if(ar.p(jp) > ar.ub(jp))
        ar.p(jp) = ar.ub(jp);
        C{jp,3} = ar.p(jp);
        set(hObject, 'Data', C);
    end
elseif(jt==5) % qLog10
    if(newval == 1)
        if(ar.qLog10(jp)==0)
            ar.qLog10(jp) = true;
            if(ar.p(jp)<=1e-10)
                ar.p(jp) = -10;
            else
                ar.p(jp) = log10(ar.p(jp));
            end
            if(ar.lb(jp)<=1e-10)
                ar.lb(jp) = -10;
            else
                ar.lb(jp) = log10(ar.lb(jp));
            end
            ar.ub(jp) = log10(ar.ub(jp));
            C{jp,2} = ar.lb(jp);
            C{jp,3} = ar.p(jp);
            C{jp,4} = ar.ub(jp);
            set(hObject, 'Data', C);
        end
    else
        if(ar.qLog10(jp)==1)
            ar.qLog10(jp) = false;
            ar.p(jp) = 10^(ar.p(jp));
            ar.lb(jp) = 10^(ar.lb(jp));
            ar.ub(jp) = 10^(ar.ub(jp));
            C{jp,2} = ar.lb(jp);
            C{jp,3} = ar.p(jp);
            C{jp,4} = ar.ub(jp);
            set(hObject, 'Data', C);
        end
    end
elseif(jt==6) % qFit
    ar.qFit(jp) = newval;
elseif(jt==7) % type
    ar.type(jp) = find(ismember(type, newval))-1;
elseif(jt==8) % mean
    newval = str2double(newval);
    ar.mean(jp) = newval;
elseif(jt==9) % std
    newval = str2double(newval);
    ar.std(jp) = newval;
end
  

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

arPlot;
