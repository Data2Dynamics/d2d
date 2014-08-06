function varargout = arNavigator(varargin)
% ARNAVIGATOR MATLAB code for arNavigator.fig
%      ARNAVIGATOR, by itself, creates a new ARNAVIGATOR or raises the existing
%      singleton*.
%
%      H = ARNAVIGATOR returns the handle to a new ARNAVIGATOR or the handle to
%      the existing singleton*.
%
%      ARNAVIGATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARNAVIGATOR.M with the given input arguments.
%
%      ARNAVIGATOR('Property','Value',...) creates a new ARNAVIGATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arNavigator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arNavigator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arNavigator

% Last Modified by GUIDE v2.5 14-Mar-2014 10:41:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arNavigator_OpeningFcn, ...
                   'gui_OutputFcn',  @arNavigator_OutputFcn, ...
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


% --- Executes just before arNavigator is made visible.
function arNavigator_OpeningFcn(hObject, eventdata, handles, varargin)
global ar
% Choose default command line output for arNavigator
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
ar.navi.index = 1;
ar.navi.C = {};
ar.navi.jm = [];
ar.navi.jplot = [];
for jm=1:length(ar.model)
    for jplot=1:length(ar.model(jm).plot)
        ar.navi.C{end+1} = ar.model(jm).plot(jplot).name;
        ar.navi.jm(end+1) = jm;
        ar.navi.jplot(end+1) = jplot;
    end
end
set(handles.popupmenu1, 'String', ar.navi.C)
refresh_all(handles)

% --- Outputs from this function are returned to the command line.
function varargout = arNavigator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
global ar;
ar.navi.index = get(hObject,'Value');
refresh_all(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global ar;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function myplot(handles, fastPlot, silent, evalfun)
arPlot(false, fastPlot, silent, evalfun, get(handles.checkbox1, 'value')==1, false);

function refresh_all(handles)
global ar;
jm = ar.navi.jm(ar.navi.index);
jplot = ar.navi.jplot(ar.navi.index);
jds = ar.model(jm).plot(jplot).dLink;
jd = jds(1);

Conds = unique(ar.model(jm).plot(jplot).condition);
qFitCond = true(size(Conds));
qFitObs = true(size(ar.model(jm).data(jd).y));

set(handles.popupmenu1, 'Value', ar.navi.index);
set(handles.uitable1, 'Data', [ar.model(jm).data(jd).y' num2cell(qFitObs)']);
set(handles.uitable2, 'Data', [Conds' num2cell(qFitCond==1)']);

pLabels = {};
for jd=jds
    ptmp = ar.model(jm).data(jd).p(ar.qDynamic(ar.model(jm).data(jd).pLink)==0);
    pLabels = union(pLabels, ptmp);
end
ar.navi.pindex = 1;
ar.navi.pLabels = pLabels;
ar.navi.jp = find(ismember(ar.pLabel, ar.navi.pLabels));
set(handles.popupmenu2, 'String', ar.navi.pLabels)

refresh_all_p(handles);


function refresh_all_p(handles)
global ar;
set(handles.popupmenu2, 'Value', ar.navi.pindex);
set(handles.edit1, 'String', ar.lb(ar.navi.jp(ar.navi.pindex)));
set(handles.edit2, 'String', ar.p(ar.navi.jp(ar.navi.pindex)));
set(handles.edit3, 'String', ar.ub(ar.navi.jp(ar.navi.pindex)));

% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
global ar

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
global ar;
jm = ar.navi.jm(ar.navi.index);
jplot = ar.navi.jplot(ar.navi.index);
jds = ar.model(jm).plot(jplot).dLink;
D1 = get(handles.uitable1, 'Data');
D2 = get(handles.uitable2, 'Data');
for jc = 1:length(jds)
    jd = jds(jc);
    qc = true;
    if(~isempty(D2))
        qc = D2{ismember(D2(:,1), ar.model(jm).plot(jplot).condition{jc}), 2};
    end
    for jo = 1:length(ar.model(jm).data(jd).y)
        ar.model(jm).data(jd).qFit(jo) = D1{jo,2} && qc;
    end
end
myplot(handles, false, true, false);

% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function uitable2_CreateFcn(hObject, eventdata, handles)
global ar

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global ar
for jm=1:length(ar.model)
    ar.model(jm).qPlotXs(:) = 0;
    ar.model(jm).qPlotYs(:) = 0;
    ar.model(jm).qPlotVs(:) = 0;
end
jm = ar.navi.jm(ar.navi.index);
jplot = ar.navi.jplot(ar.navi.index);
ar.model(jm).qPlotYs(jplot) = 1;
myplot(handles, false, true, false);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ar
ar.navi.index = ar.navi.index - 1;
if(ar.navi.index<1)
    ar.navi.index = 1;
end
refresh_all(handles)
for jm=1:length(ar.model)
    ar.model(jm).qPlotXs(:) = 0;
    ar.model(jm).qPlotYs(:) = 0;
    ar.model(jm).qPlotVs(:) = 0;
end
jm = ar.navi.jm(ar.navi.index);
jplot = ar.navi.jplot(ar.navi.index);
ar.model(jm).qPlotYs(jplot) = 1;
myplot(handles, false, true, false);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ar
ar.navi.index = ar.navi.index + 1;
if(ar.navi.index>length(ar.navi.C))
    ar.navi.index = length(ar.navi.C);
end
refresh_all(handles)
for jm=1:length(ar.model)
    ar.model(jm).qPlotXs(:) = 0;
    ar.model(jm).qPlotYs(:) = 0;
    ar.model(jm).qPlotVs(:) = 0;
end
jm = ar.navi.jm(ar.navi.index);
jplot = ar.navi.jplot(ar.navi.index);
ar.model(jm).qPlotYs(jplot) = 1;
myplot(handles, false, true, false);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
global ar;
ar.navi.pindex = get(hObject,'Value');
refresh_all_p(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global ar
ar.navi.pindex = ar.navi.pindex - 1;
if(ar.navi.pindex<1)
    ar.navi.pindex = 1;
end
refresh_all_p(handles)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global ar
ar.navi.pindex = ar.navi.pindex + 1;
if(ar.navi.pindex>length(ar.navi.pLabels))
    ar.navi.pindex = length(ar.navi.pLabels);
end
refresh_all_p(handles)


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
global ar;
ip = ar.navi.jp(ar.navi.pindex);
ar.p(ip) = ar.p(ip) - 0.1;
if(ar.p(ip) < ar.lb(ip))
    ar.p(ip) = ar.lb(ip);
end
if(ar.p(ip) > ar.ub(ip))
    ar.p(ip) = ar.ub(ip);
end
% arChi2(false)
myplot(handles, true, true, true);
refresh_all_p(handles);


function edit2_Callback(hObject, eventdata, handles)
global ar;
ip = ar.navi.jp(ar.navi.pindex);
newval = str2double(get(hObject,'String'));
if(~isnan(newval))
    if(newval < ar.lb(ip))
        newval = ar.lb(ip);
    end
    if(newval > ar.ub(ip))
        newval = ar.ub(ip);
    end
    ar.p(ip) = newval;
end
set(hObject, 'String', ar.p(ip));
% arChi2(false)
myplot(handles, true, true, true);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
global ar;
ip = ar.navi.jp(ar.navi.pindex);
ar.p(ip) = ar.p(ip) + 0.1;
if(ar.p(ip) < ar.lb(ip))
    ar.p(ip) = ar.lb(ip);
end
if(ar.p(ip) > ar.ub(ip))
    ar.p(ip) = ar.ub(ip);
end
% arChi2(false)
myplot(handles, true, true, true);
refresh_all_p(handles);

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ar;
arFitSingle(ar.navi.jp(ar.navi.pindex), true);
%arChi2(false)
myplot(handles, true, true, true);
refresh_all_p(handles);
