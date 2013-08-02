function varargout = arTuner(varargin)
% ARTUNER M-file for arTuner.fig
% GUI to assign parameter values

% Last Modified by GUIDE v2.5 25-Jun-2013 11:48:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arTuner_OpeningFcn, ...
                   'gui_OutputFcn',  @arTuner_OutputFcn, ...
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

% --- Executes just before arTuner is made visible.
function arTuner_OpeningFcn(hObject, eventdata, handles, varargin)
global ar;
% Choose default command line output for arTuner
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
ar.tuner.index = 1;

% --- Outputs from this function are returned to the command line.
function varargout = arTuner_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
global ar;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', ar.pLabel)

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
global ar;
ar.tuner.index = get(hObject,'Value');
refresh_all(handles)

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global ar;
ar.p(ar.tuner.index) = ar.p(ar.tuner.index) + 0.1;
if(ar.p(ar.tuner.index) < ar.lb(ar.tuner.index))
    ar.p(ar.tuner.index) = ar.lb(ar.tuner.index);
end
if(ar.p(ar.tuner.index) > ar.ub(ar.tuner.index))
    ar.p(ar.tuner.index) = ar.ub(ar.tuner.index);
end
% arChi2(false)
arPlot(false, true, true)
refresh_textbox(handles)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global ar;
ar.p(ar.tuner.index) = ar.p(ar.tuner.index) - 0.1;
if(ar.p(ar.tuner.index) < ar.lb(ar.tuner.index))
    ar.p(ar.tuner.index) = ar.lb(ar.tuner.index);
end
if(ar.p(ar.tuner.index) > ar.ub(ar.tuner.index))
    ar.p(ar.tuner.index) = ar.ub(ar.tuner.index);
end
% arChi2(false)
arPlot(false, true, true)
refresh_textbox(handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
global ar;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', ar.p(1));

function edit1_Callback(hObject, eventdata, handles)
global ar;
newval = str2double(get(hObject,'String'));
if(~isnan(newval))
    if(newval < ar.lb(ar.tuner.index))
        newval = ar.lb(ar.tuner.index);
    end
    if(newval > ar.ub(ar.tuner.index))
        newval = ar.ub(ar.tuner.index);
    end
    ar.p(ar.tuner.index) = newval;
end
set(hObject, 'String', ar.p(ar.tuner.index));
% arChi2(false)
arPlot(false, true)

% --- Executes during object creation, after setting all properties.
function checkbox1_CreateFcn(hObject, eventdata, handles)
global ar;
if(ar.qLog10(1)==1)
    set(hObject, 'Value', 1);
else
    set(hObject, 'Value', 0);
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
global ar;
if(get(hObject,'Value') == 1)
    if(ar.qLog10(ar.tuner.index)==0)
        ar.qLog10(ar.tuner.index) = true;
        if(ar.p(ar.tuner.index)<=1e-10)
            ar.p(ar.tuner.index) = -10;
        else
            ar.p(ar.tuner.index) = log10(ar.p(ar.tuner.index));
        end
        if(ar.lb(ar.tuner.index)<=1e-10)
            ar.lb(ar.tuner.index) = -10;
        else
            ar.lb(ar.tuner.index) = log10(ar.lb(ar.tuner.index));
        end
		ar.ub(ar.tuner.index) = log10(ar.ub(ar.tuner.index));
    end
else
    if(ar.qLog10(ar.tuner.index)==1)
        ar.qLog10(ar.tuner.index) = false;
        ar.p(ar.tuner.index) = 10^(ar.p(ar.tuner.index));
		ar.lb(ar.tuner.index) = 10^(ar.lb(ar.tuner.index));
		ar.ub(ar.tuner.index) = 10^(ar.ub(ar.tuner.index));
    end
end
refresh_all(handles)

% --- Executes during object creation, after setting all properties.
function radiobutton4_CreateFcn(hObject, eventdata, handles)
global ar;
if(ar.qFit(1)==1)
    set(hObject, 'Value', 1);
end

% --- Executes during object creation, after setting all properties.
function radiobutton5_CreateFcn(hObject, eventdata, handles)
global ar;
if(ar.qFit(1)==0)
    set(hObject, 'Value', 1);
end

% --- Executes during object creation, after setting all properties.
function radiobutton6_CreateFcn(hObject, eventdata, handles)
global ar;
if(ar.qFit(1)==2)
    set(hObject, 'Value', 1);
end

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
global ar;
if(get(hObject,'Value') == 1)
    ar.qFit(ar.tuner.index) = 1;
end

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
global ar;
if(get(hObject,'Value') == 1)
    ar.qFit(ar.tuner.index) = 0;
end

% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
global ar;
if(get(hObject,'Value') == 1)
    ar.qFit(ar.tuner.index) = 2;
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
arFit(true);
arChi2(false)
arPlot(false, true)
refresh_textbox(handles)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
arFitObs(true);
arChi2(false)
arPlot(false, true)
refresh_textbox(handles)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
arFitDyn(true);
arChi2(false)
arPlot(false, true)
refresh_textbox(handles)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
global ar;
arFitInit(true);
arChi2(false)
arPlot(false, true)
refresh_textbox(handles)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global ar;
arFitSingle(ar.tuner.index, true);
arChi2(false)
arPlot(false, true)
refresh_textbox(handles)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
global ar;
ar.tuner.index = ar.tuner.index - 1;
if(ar.tuner.index <= 0)
    ar.tuner.index = 1;
else
    refresh_all(handles)
end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
global ar;
ar.tuner.index = ar.tuner.index + 1;
if(ar.tuner.index > length(ar.p))
    ar.tuner.index = length(ar.p);
else
    refresh_all(handles)
end



function edit3_Callback(hObject, eventdata, handles)
global ar;
newval = str2double(get(hObject,'String'));
if(~isnan(newval))
    ar.ub(ar.tuner.index) = newval;
    if(ar.p(ar.tuner.index) > ar.ub(ar.tuner.index))
        ar.p(ar.tuner.index) = ar.ub(ar.tuner.index);
    end
end
set(hObject, 'String', ar.ub(ar.tuner.index));
refresh_textbox(handles)

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
global ar;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', ar.ub(1));


function edit4_Callback(hObject, eventdata, handles)
global ar;
newval = str2double(get(hObject,'String'));
if(~isnan(newval))
    ar.lb(ar.tuner.index) = newval;
    if(ar.p(ar.tuner.index) < ar.lb(ar.tuner.index))
        ar.p(ar.tuner.index) = ar.lb(ar.tuner.index);
    end
end
set(hObject, 'String', ar.lb(ar.tuner.index));
refresh_textbox(handles)


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
global ar;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', ar.lb(1));

%% Helper Functions
function refresh_textbox(handles)
global ar;
set(handles.edit1, 'String', ar.p(ar.tuner.index));

function refresh_all(handles)
global ar;
set(handles.popupmenu1, 'Value', ar.tuner.index);
set(handles.edit1, 'String', ar.p(ar.tuner.index));
set(handles.edit3, 'String', ar.ub(ar.tuner.index));
set(handles.edit4, 'String', ar.lb(ar.tuner.index));
if(ar.qLog10(ar.tuner.index)==1)
    set(handles.checkbox1, 'Value', 1);
else
    set(handles.checkbox1, 'Value', 0);
end
if(ar.qFit(ar.tuner.index)==1)
    set(handles.radiobutton4, 'Value', 1);
elseif(ar.qFit(ar.tuner.index)==0)
    set(handles.radiobutton5, 'Value', 1);
elseif(ar.qFit(ar.tuner.index)==2)
    set(handles.radiobutton6, 'Value', 1);
end

