function varargout = arPlotter(varargin)
% ARPLOTTER MATLAB code for arPlotter.fig
%      ARPLOTTER, by itself, creates a new ARPLOTTER or raises the existing
%      singleton*.
%
%      H = ARPLOTTER returns the handle to a new ARPLOTTER or the handle to
%      the existing singleton*.
%
%      ARPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARPLOTTER.M with the given input arguments.
%
%      ARPLOTTER('Property','Value',...) creates a new ARPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arPlotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arPlotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arPlotter

% Last Modified by GUIDE v2.5 19-Aug-2013 07:52:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arPlotter_OpeningFcn, ...
                   'gui_OutputFcn',  @arPlotter_OutputFcn, ...
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


% --- Executes just before arPlotter is made visible.
function arPlotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to arPlotter (see VARARGIN)

% Choose default command line output for arPlotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes arPlotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = arPlotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)

global ar

% column name
C = get(hObject, 'ColumnName');
if(ar.config.fiterrors == 1)
    C{2} = '-2*log(L)';
else
    C{2} = 'chi^2';
end
set(hObject,'ColumnName',C);

% table data
% get(hObject)
ar.plotter.C = {};
ar.plotter.jm = [];
ar.plotter.jplot = [];
llhreltodata = [];
for jm=1:length(ar.model)
    for jplot=1:length(ar.model(jm).plot)
        chi2 = 0;
        ndata = 0;
        logfitting = 0;
        logplotting = 0;
        qfit = 0;
        if(isfield(ar.model(jm),'data'))
            for jd = ar.model(jm).plot(jplot).dLink
                chi2 = chi2 + sum(ar.model(jm).data(jd).chi2);
                if(ar.config.fiterrors == 1)
                    chi2 = chi2 + sum(ar.model(jm).data(jd).chi2err);
                end
                ndata = ndata + sum(ar.model(jm).data(jd).ndata);
                logfitting = logfitting + sum(ar.model(jm).data(jd).logfitting);
                logplotting = logplotting + sum(ar.model(jm).data(jd).logplotting);
                qfit = qfit + sum(ar.model(jm).data(jd).qFit==0);
            end
            
            if(ar.config.fiterrors == 1)
                chi2 = 2*ndata*log(sqrt(2*pi)) + chi2;
            end
        end
        llhreltodata(end+1) = chi2/ndata; %#ok<AGROW>
        ar.plotter.C{end+1,1} = ar.model(jm).plot(jplot).name;
        ar.plotter.C{end,2} = chi2;
        ar.plotter.C{end,3} = ndata;
        ar.plotter.C{end,4} = length(ar.model(jm).plot(jplot).dLink);
        ar.plotter.C{end,5} = ar.model(jm).qPlotYs(jplot);
        ar.plotter.C{end,6} = ar.model(jm).qPlotXs(jplot);
        ar.plotter.C{end,7} = ar.model(jm).qPlotVs(jplot);
        ar.plotter.C{end,8} = logplotting>0;
        ar.plotter.C{end,9} = logfitting>0;
        ar.plotter.C{end,10} = qfit==0;
        
        ar.plotter.jm(end+1) = jm;
        ar.plotter.jplot(end+1) = jplot;
    end
end
if(isfield(ar.config,'plotter_sorted') && ar.config.plotter_sorted)
    [~, isorted] = sort(llhreltodata,2,'descend');
    ar.plotter.C = ar.plotter.C(isorted,:);
    ar.plotter.jm = ar.plotter.jm(isorted);
    ar.plotter.jplot = ar.plotter.jplot(isorted);
end

set(hObject, 'Data', ar.plotter.C);


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

global ar

jm = ar.plotter.jm(eventdata.Indices(1));
jplot = ar.plotter.jplot(eventdata.Indices(1));

if(eventdata.Indices(2)==5) % PlotY
    ar.model(jm).qPlotYs(jplot) = eventdata.EditData;
elseif(eventdata.Indices(2)==6) % PlotX
    ar.model(jm).qPlotXs(jplot) = eventdata.EditData;
elseif(eventdata.Indices(2)==7) % PlotV
    ar.model(jm).qPlotVs(jplot) = eventdata.EditData;
elseif(eventdata.Indices(2)==8) % PlotLog
    for jd = ar.model(jm).plot(jplot).dLink
        ar.model(jm).data(jd).logplotting(:) = eventdata.EditData;
    end
elseif(eventdata.Indices(2)==9) % FitLog
    ar.plotter.C{eventdata.Indices(1), 8} = eventdata.EditData;
    ar.plotter.C{eventdata.Indices(1), 9} = eventdata.EditData;
    set(hObject, 'Data', ar.plotter.C);
    for jd = ar.model(jm).plot(jplot).dLink
        if(~eventdata.EditData)
            ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==1) = ...
                10.^ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==1);
        else
            ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==0) = ...
                log10(ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==0));
        end
        ar.model(jm).data(jd).logfitting(:) = eventdata.EditData;
    end
elseif(eventdata.Indices(2)==10) % Fit
    for jd = ar.model(jm).plot(jplot).dLink
        ar.model(jm).data(jd).qFit(:) = eventdata.EditData;
    end
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
arPlot(false, false, true);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all
