function nomadm_functions(choice,param);
%NOMADM_FUNCTIONS  Run NOMADm GUI functions.
%
%   Syntax:
%      nomadm_functions(CHOICE)
%      nomadm_functions(CHOICE,PARAM)
%
%   Description:
%      NOMADM_FUNCTIONS executes one of several NOMADm
%      GUI functions.  CHOICE is a string that describes
%      which GUI function is to be performed.  Most of
%      the GUI functions are specified in callback
%      functions associated with GUI objects, such as
%      menu choices, pushbuttons, etc., as defined by
%      the calling function NOMADm.  The optional
%      variable PARAM is used to pass additional data
%      needed by an individual GUI function.
%
%   Note:  This function should only be called by NOMADM 
%   and not run independently.
%
%   See also NOMADM, MADS

%**************************************************************************
% Copyright (c) 2001-2004 by Mark A. Abramson
%
% This file is part of the NOMADm software package.
%
% NOMADm is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% NOMADm is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% for more details.
%
% You should have received a copy of the GNU General Public License along
% with NOMADm; if not, write to the Free Software Foundation, Inc., 
% 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% -------------------------------------------------------------------------
% Author information:
%   Mark A. Abramson, LtCol, USAF, PhD
%   Air Force Institute of Technology
%   Department of Mathematics and Statistics
%   2950 Hobson Way
%   Wright-Patterson AFB, OH 45433
%   (937) 255-3636 x4524
%   Mark.Abramson@afit.edu
%**************************************************************************


%**************************************************************************
% NOMADm_Functions:  Run GUI functions associated with NOMADm
% -------------------------------------------------------------------------
% Called by: NOMADm
% Calls:     runMADS,       loadMADS,    updateScreen, 
%            searchOptions, daceOptions, viewResults
% VARIABLES:
%  choice     = string defining user-selected menu choice
%  param      = additional parameter passed with choice
%  GUI        = structure containing all GUI handles and variables
%  onoff      = cell array with two strings "on" and "off"
%  filterSpec = cell array of file types for file input dialog boxes
%  RunSet     = global variable used to store MADS output data
%
%  All other variable specified below before menu choice
%**************************************************************************

global GUI
onoff = {'on','off'};
filterSpec = {'*.m',        'Matlab M-files (*.m)'; ...
              '*.f; *.for', 'Fortran files (*.f,*.for)'; ...
              '*.c; *.C',   'C/C++ files (*.c,*.C)'};
if isappdata(0,'RUNSET'), RunSet = getappdata(0,'RUNSET'); end

% Perform only the operations indicated by the "choice" variable
switch choice

   %***********************************************************************
   % ChangeProblem: Choose an optimization problem
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI               = structure of all GUI handles and variables
   %    .Path           =   path of current optimization problem
   %    .ProblemExt     =   filename extension of optimization problem 
   %    .ProblemName    =   name of current optimization problem
   %    .RunStatus      =   handle for text bar at bottom of figure window
   %    .RunMenuXXX     =   handles for Run menu items
   %    .ProblemMenuXXX =   handles for Problem menu items
   %    .MADSMenuXXX    =   handles for MADS menu items
   %    .SessionMenuXXX =   handles for Session menu items
   %    .CacheMenuXXX   =   handles for CacheMenu items
   %    .FileExt        =   filename suffixes
   %      .O            =     Omega file suffix
   %      .I            =     initial points file suffix
   %      .N            =     discrete neighbors file suffix
   %      .P            =     user parameter file suffix
   %      .C            =     Cache file suffix
   %      .S            =     Session file suffix
   %    .NoGrad         =   flag indicating initial availability of gradients
   %  PFile             = filename of current optimization problem
   %  PPath             = pathname of current optimization problem
   %  ProblemName       = name of current optimization problem
   %  isConstrained     = flags problem as having linear constraints
   %  hasInitGuess      = flags problem as having an initial guess
   %  hasUserParam      = flags problem as having user parameters
   %  isMVP             = flags problem as being an MVP
   %  existCache        = flags problem as having an existing Cache file
   %  existSession      = flags problem as having an existing Session file
   %  PLabels           = strings of Poll strategy labels
   %***********************************************************************
case 'ChangeProblem'
   [PFile,PPath] = uigetfile(filterSpec, ...
                             'Select an Optimization Problem file');
   if (PFile)
      nomadm_functions('Clear');
      nomadm_functions('ResetSession');
      rmpath(GUI.Path);
      GUI.Path = PPath;
      addpath(GUI.Path);
      [ProblemName, GUI.ProblemExt] = strtok(PFile,'.');
      set(GUI.ProblemName,'String',ProblemName);
      set(GUI.RunStatus,  'String','No runs performed');
      set(GUI.RunMenuExec,        'Enable','on');
      set(GUI.RunMenuOneIteration,'Enable','on');
      set(GUI.RunMenuExecFeasible,'Enable','on');
      
      % Allow editing only of files that exist
      isConstrained = exist([ProblemName,GUI.FileExt.O,'.m']);
      hasInitGuess  = exist([ProblemName,GUI.FileExt.I,'.m']);
      isMVP         = exist([ProblemName,GUI.FileExt.N,'.m']);
      hasUserParam  = exist([ProblemName,GUI.FileExt.P,'.m']);
      existCache    = exist([ProblemName,GUI.FileExt.C]);
      existSession  = exist([ProblemName,GUI.FileExt.S]);
      set(GUI.ProblemMenuEdit(1),     'Visible', 'on');
      set(GUI.ProblemMenuEdit(2),     'Visible', onoff{1+~hasInitGuess});
      set(GUI.ProblemMenuEdit(3),     'Visible', onoff{1+~isConstrained});
      set(GUI.ProblemMenuEdit(4),     'Visible', onoff{1+~isMVP});
      set(GUI.ProblemMenuEdit(5),     'Visible', onoff{1+~hasUserParam});
      set(GUI.MADSMenuPollOrder(end), 'Enable',  onoff{1+~hasUserParam});
      set(GUI.SessionMenuLoad,        'Enable',  onoff{1+~existSession});
      set(GUI.SessionMenuSave,        'Enable',  'on');
      set(GUI.SessionMenuDelete,      'Enable',  onoff{1+~existSession});
      set(GUI.CacheMenuLoad,          'Enable',  onoff{1+~existCache});
      set(GUI.CacheMenuCount,         'Enable',  onoff{1+~existCache});
      set(GUI.CacheMenuDelete,        'Enable',  onoff{1+~existCache});      
      
      % Allow gradient options if functions file returns enough arguments
      try
         GUI.NoGrad = (abs(nargout(ProblemName)) <= 2);
      catch
         GUI.NoGrad = 0;
      end
      PLabels = char(GUI.Labels.PollStrategy);
      set(GUI.MADSMenuPollStrategy(find(PLabels(:,1) == 'G')), ...
                                  'Enable',onoff{1+GUI.NoGrad});
   end

   %***********************************************************************
   % EditFile: Edit an optimization problem file
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  k              = Problem Menu position number of selected file
   %  Ext            = filename extension of file to be edited
   %  FileName       = full file name of file to be edited
   %  GUI            = structure of all GUI handles and variables
   %    .Path        =   path of current optimization problem
   %    .ProblemExt  =   filename extension of optimization problem 
   %    .ProblemName =   name of current optimization problem
   %    .Types.File  = strings of file suffixes
   %***********************************************************************
case 'EditFile'
   k = get(gcbo,'Position') - 1;
   if (k == 1)
      Ext = GUI.ProblemExt;
   else
      Ext = '.m';
   end
   FileName = fullfile(GUI.Path, ...
                  [get(GUI.ProblemName,'String'),GUI.Types.File{k},Ext]);
   edit(FileName);
   clear(FileName);

   %***********************************************************************
   % Scaling: Set the variable scaling of an optimization problem
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI.ScaleMenu = handle for GUI Scaling submenu of MADS menu
   %***********************************************************************
case 'Scaling'
   set(GUI.ScaleMenu,'Checked','off');
   set(GUI.ScaleMenu(get(gcbo,'Position')),'Checked','on');   

   %***********************************************************************
   % LoadSession: Load session options from previously saved file
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     updateScreen
   % VARIABLES:
   %  SessionFile     = full path and name of session file
   %  GUI.Path        = path of current optimization problem
   %  GUI.ProblemName = name of current optimization problem
   %  GUI.FileExt.S   = filename extension of session file
   %  GUI.Choice      = current user choices
   %  GUI.Options     = current option settings
   %  Session.Choice  = previously saved user choices 
   %  Session.Options = previously saved option settings
   %***********************************************************************
case 'LoadSession'
   SessionFile = fullfile(GUI.Path, ...
                 [get(GUI.ProblemName,'String'), GUI.FileExt.S]);
   if exist(SessionFile,'file')
      load(SessionFile);
      GUI.Choice  = Session.Choice;
      GUI.Options = Session.Options;
      GUI = updateScreen('all', GUI.Choice, GUI.Options, GUI);
   end

   %***********************************************************************
   % ResetSession: Reset MADS Parameters to Default values
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     updateScreen
   % VARIABLES:
   %  GUI.Choice           = current user choices
   %  GUI.Options          = current option settings
   %  GUI.Defaults.Choice  = default user choice settings
   %  GUI.Defaults.Options = default options settings
   %***********************************************************************
case 'ResetSession'
   GUI.Choice  = GUI.Defaults.Choice;
   GUI.Options = GUI.Defaults.Options;
   GUI = updateScreen('all', GUI.Choice, GUI.Options, GUI);

   %***********************************************************************
   % EditSearch: Edit Search Parameters
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     searchOptions
   % VARIABLES:
   %  NewAns                    = answers to an input dialog box
   %  GUI.Options.nSearches     = number of Search types to be used
   %  GUI.Options.SurrOptimizer = string ID for surrogate optimizer
   %  GUI.SearchType            = handles for displayed Search types
   %  nSearches                 = number of Search types selected
   %***********************************************************************
case 'EditSearch'
   NewAns = inputdlg({'Number of Search types','Surrogate Optimizer'}, ...
                      'Search Information', 1, ...
                      {int2str(GUI.Options.nSearches), ...
                       GUI.Options.SurrOptimizer});
   if isempty(NewAns), return, end
   nSearches = str2num(NewAns{1});
   GUI.Options.SurrOptimizer = NewAns{2};
   if isempty(nSearches) | length(nSearches) > 1, return, end
   if (nSearches == 0)
      set(GUI.SearchType(1:end), 'String', 'None');
   else
      GUI.Options.nSearches = nSearches;
      searchOptions(nSearches,GUI.Options.SurrOptimizer);
   end

   %***********************************************************************
   % LoadSearchOptions: Load selected Search parameters
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     updateScreen
   % VARIABLES:
   %  k               = Search number
   %  GUI             = structure of all GUI handles and variables
   %    .Choice       =   structure of user choices
   %      .search(k)  =     user Search choices
   %    .Options      =   structure of MADS parameters settings
   %      .Search(k)  =     structure of user-selected k-th Search
   %      .nSearches  = number of Search types to be used
   %    .Types.Search = list of possible Search types
   %    .figSearch    = handle for the Search figure window
   %    .SearchScreen = handles for Search figure window objects
   %    .NoGrad       =   flag indicating initial availability of gradients
   %  Search          = temporary storage of GUI.Options.Search(k)
   %    .type         =   string identifying Search type
   %    .label        =   long text label for Search type
   %    .nIter        =   number of iterations to perform Search
   %    .nPoints      =   number of Search points
   %    .file         =   optional user file defining Search
   %***********************************************************************
case 'LoadSearchOptions'
   loadSearch = get(gcbo,'UserData');
   if loadSearch
      for k = 1:GUI.Options.nSearches
         GUI.Choice.search(k) = get(GUI.SearchScreen(k).type, 'Value');
         Search.type = GUI.Types.Search{GUI.Choice.search(k)};
         if (GUI.NoGrad & strcmp(Search.type, 'GPollI'))
            uiwait(msgbox('No derivatives available for this problem', ...
                         ['Error in Search Type #', int2str(k)], ...
                          'error','modal'));
            return
         end
         Search.label   = GUI.Labels.Search{GUI.Choice.search(k)};
         Search.nIter   = str2num( get(GUI.SearchScreen(k).nIter,  'String'));
         Search.nPoints = str2num( get(GUI.SearchScreen(k).nPoints,'String'));
         Search.file    = fullfile(get(GUI.SearchScreen(k).file, 'UserData'), ...   
                                get(GUI.SearchScreen(k).file, 'String'));
         GUI.Options.Search(k) = Search;
      end
      GUI = updateScreen('search', GUI.Choice,GUI.Options,GUI);
   end
   close(GUI.figSearch);

   %***********************************************************************
   % LoadDACEOptions: Select a DACE surrogate option
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  k                  = Search number
   %  dace               = structure of DACE parameters
   %    .Reg             =   handle for regression function
   %    .Corr            =   handle for correlation function
   %    .Theta           =   initial guess for correlation parameters
   %    .Lower           =   lower bounds for correlation parameters
   %    .Upper           =   upper bounds for correlation parameters
   %    .Isotropic       =   flag for isotropic correlation parameters
   %  GUI                = structure of all GUI handles and variables
   %    .Types           =   lists of possible types
   %      .daceReg       =     list of possible DACE regression functions
   %      .daceCorr      =     list of possible DACE correlation functions
   %    .SearchScreen    =   handles for Search figure window objects
   %      .file          =     optional user file defining Search
   %    .daceOptionsXXX  =   handles for DACE figure window objects
   %    .Options.dace(k) =   structure of DACE parameters for Search k
   %    .figDACE         =   handle for DACE figure window
   %***********************************************************************
case 'LoadDACEOptions'
   k = get(gcbo, 'UserData');
   if k
      dace.Reg      = GUI.Types.daceReg{ get(GUI.daceOptionsReg, 'Value')};
      dace.Corr     = GUI.Types.daceCorr{get(GUI.daceOptionsCorr,'Value')};
      dace.Theta    = str2num(get(GUI.daceOptionsTheta,'String'));
      dace.Lower    = str2num(get(GUI.daceOptionsLower,'String'));
      dace.Upper    = str2num(get(GUI.daceOptionsUpper,'String'));
      dace.Isotropic = get(GUI.daceOptionsIsotropic,'Value');
      set(GUI.SearchScreen(k).file, 'UserData', '');
      set(GUI.SearchScreen(k).file, 'String', dace.Reg);
      GUI.Options.dace(k) = dace;      
   end
   close(GUI.figDACE);

   %***********************************************************************
   % LoadUserSearchFile: Load a user-specified Search file
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     daceOptions
   % VARIABLES:
   %  k                     = Search number
   %  GUI.Types.Search      = list of possible Search types
   %  GUI.SearchScreen.file = handle for Search figure window user file
   %  SName, SPath          = name and path of user Search file
   %  filterSpec            = file types for use in input dialog boxes
   %***********************************************************************
case 'LoadUserSearchFile'
   k = get(gcbo, 'UserData');
   switch GUI.Types.Search{get(gcbo, 'Value')}
   case {'Custom','CustomS'}
      [SName,SPath] = uigetfile(filterSpec, ...
                               ['Choose File for Search #', int2str(k)]);
      if (SName)
         set(GUI.SearchScreen(k).file, 'UserData', SPath);
         set(GUI.SearchScreen(k).file, 'String',   SName);
      end
   case {'DACE'}
      daceOptions(k);
   otherwise
      set(GUI.SearchScreen(k).file, 'String', '');
   end

   %***********************************************************************
   % Edit Parameters: Edit parameters
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  param   = User data storage
   %  Title   = title of input dialog box
   %  Labels  = text labels in input dialog box
   %  Handles = handles to appropriate GUI figure objects
   %  DefAns  = default answers that appear in input dialog box
   %  NewAns  = new answers that appear in input dialog box
   %***********************************************************************
case 'EditParameters'
   param   = get(gcbo,'UserData');
   [Title, Labels, Handles] = deal(param{:});
   DefAns  = get(Handles, 'String');
   NewAns  = inputdlg(Labels,Title,ones(length(DefAns),1)*[1,45],DefAns);
   if (~isempty(NewAns))
      set(Handles, {'String'}, NewAns);
   end

   %***********************************************************************
   % SelectPollXXX: Select a Poll Strategy, Poll Center, or Poll Order
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI.Choice      = user choices
   %    .pollStrategy =   selected Poll strategy
   %    .pollCenter   =   selected Poll center
   %    .pollOrder    =   selected Poll order strategy
   %  GUI.pollXXX     = GUI figure window objects for GUI.Choice items
   %***********************************************************************
case 'SelectPollStrategy'
   GUI.Choice.pollStrategy = get(gcbo,'Position');
   set(GUI.PollStrategy,'String', get(gcbo,'Label'));
case 'SelectPollCenter'
   GUI.Choice.pollCenter = get(gcbo,'Position');
   set(GUI.PollCenter,'String', get(gcbo,'Label'));
case 'SelectPollOrder'
   GUI.Choice.pollOrder = get(gcbo,'Position');
   set(GUI.PollOrder,'String', get(gcbo,'Label'));

   %***********************************************************************
   % Toggle: Toggle the check mark of a menu item
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   %***********************************************************************
case 'Toggle'
   umtoggle(gcbo);

   %***********************************************************************
   % SaveCache: Save MADS run results in a .mat file for later retrieval
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI             = structure of all GUI handles and variables
   %    .fig          =   handle for the GUI figure window
   %    .RunStatus    =   handle for the GUI figure window status bar
   %    .RunCount     =   MADS run number
   %    .ProblemName  =   name of current optimization problem
   %    .Path         =   path of current optimization problem
   %    .CacheMenuXXX =   handles for Cache menu items
   %    .FileExt.C    =   default Cache filename suffix 
   %  CName           = name of Cache file
   %  RunSet(1).Cache = structure of MADS run data
   %  Cache           = storage of Cache for use by MADS
   %***********************************************************************
case 'SaveCache'
   setptr(GUI.fig,'watch');
   set(GUI.RunStatus, 'String', ...
       ['Run # ',int2str(GUI.RunCount),' Cache being saved']);
   CName = [get(GUI.ProblemName, 'String'), GUI.FileExt.C];
   Cache = RunSet(1).Cache;
   save([GUI.Path, CName],'Cache');
   set([GUI.CacheMenuLoad; GUI.CacheMenuCount;GUI.CacheMenuDelete], ...
       'Enable','on');
   set(GUI.RunStatus, 'String', ...
       ['Run # ',int2str(GUI.RunCount),' Cache saved in ',CName]);
   setptr(GUI.fig,'arrow');

   %***********************************************************************
   % SaveSession: Save session options to file for future retrieval
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     loadMADS
   % VARIABLES:
   %  SessionFile       = full path and name of session file
   %  GUI               = structure of all GUI handles and variables
   %    .Path           =   path of current optimization problem
   %    .ProblemName    =   name of current optimization problem
   %    .FileExt.S      =   filename extension of session file
   %    .Choice         =   current user choices
   %    .Options        =   current option settings
   %    .SessionMenuXXX =   handles for Session menu items
   %    .OptionsMenuXXX =   handles for Options menu items
   %    .RunStatus      =   handle for figure window status bar
   %  Session           = previously saved options and parameters
   %    .Problem        =   optimization problem data
   %    .Choice         =   user choices 
   %    .Options        =   option settings
   %      .Term         =     termination criteria
   %        .iter       =       number of iterations
   %        .func       =       number of function evaluations
   %        .time       =       CPU time
   %        .fails      =       number of consecutive Poll failures
   %      .TermRel      =     termination relative to initial mesh size
   %***********************************************************************
case 'SaveSession'
   SessionFile = [get(GUI.ProblemName,'String'), GUI.FileExt.S];
   [Session.Problem,Session.Options] = loadMADS(0,GUI);
   Session.Options.Term.iter  = get(GUI.TermIter,  'String');
   Session.Options.Term.func  = get(GUI.TermFunc,  'String');
   Session.Options.Term.time  = get(GUI.TermTime,  'String');
   Session.Options.Term.fails = get(GUI.TermFails, 'String');
   umtoggle(GUI.OptionsMenuTermRel);
   Session.Options.TermRel   = umtoggle(GUI.OptionsMenuTermRel);
   Session.Choice = GUI.Choice;
   save([GUI.Path, SessionFile],'Session');
   set(GUI.RunStatus,'String', ...
       ['Session Options saved to file, ' SessionFile]);
   set([GUI.SessionMenuLoad; GUI.SessionMenuDelete],'Enable','on');

   %***********************************************************************
   % DeleteFile: Delete a Session or Cache file
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  param          = user data that identifies file to be deleted
   %  FileID         = string identifing file to be deleted
   %  FileExt        = filename suffix of file to be deleted
   %  handles        = handles of menu choices to be disabled
   %  file           = full path and name of file to be deleted
   %  deleteFile     = logical for deleting file
   %  GUI            =   structure of all GUI handles and variables
   %    .Path        =   path of current optimization problem
   %    .ProblemName =   name of current optimization problem
   %    .RunStatus   =   handle for figure window status bar
   %***********************************************************************
case 'DeleteFile'
   param = get(gcbo,'UserData');
   [FileID,FileExt,handles] = deal(param{:});
   file  = fullfile(GUI.Path,[get(GUI.ProblemName,'String'), FileExt]);
   deleteFile = questdlg(['Are you sure you want to delete ',file,'?'], ...
                         ['Delete ',FileID,'?'],'No');
   if strcmp(deleteFile, 'Yes')
      delete(file);
      set(GUI.RunStatus,'String', [FileID,', ',file,', has been deleted']);
      set(handles,'Enable','off');
   end
   
   %***********************************************************************
   % Clear: Clear figure window and reset all the appropriate variables
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI               = structure of all GUI handles and variables
   %    .fig            =   handle for GUI figure window
   %    .StopRun        =   handle for Stop Run pushbutton
   %    .ResumeRun      =   handle for Resume Run pushbutton
   %    .axesHistory    =   handle for History plot
   %    .axesFilter     =   handle for Filter plot
   %    .ResultsMenuXXX =   handles for Results menu items
   %    .RunMenuXXX     =   handles for Run menu items
   %    .CacheMenuXXX   =   handles for Cache menu items
   %    .RunStatus      =   handle for figure window status bar
   %    .RunMax         =   maximum nmber of MADS runs
   %    .RunCount       =   MADS run counter
   %  RunSet            = structure of MADS run data
   %***********************************************************************
case 'Clear'
   cla;
   setptr(GUI.fig,'arrow');
   set([GUI.axesHistory; get(GUI.axesHistory,'Children')],'Visible','off');
   set([GUI.axesFilter;  get(GUI.axesFilter, 'Children')],'Visible','off');
   set([GUI.StopRun; GUI.ResumeRun],                      'Visible','off');      
   set([GUI.ResultsMenuView(1:GUI.RunMax)'],              'Visible','off');
   set([GUI.ResultsMenu],                                 'Visible','off');
   set([GUI.RunMenuClear;  GUI.CacheMenuSave],            'Enable', 'off');
   set([GUI.RunMenuResume; GUI.RunMenuRestart],           'Enable', 'off');      
   set(GUI.axesHistory, 'NextPlot','replace');
   set(GUI.RunStatus,'String','All runs cleared');
   GUI.RunCount = 0;
   RunSet = [];
   clear functions;
   if isappdata(0,'RUNSET'), rmappdata(0,'RUNSET'); end

   %***********************************************************************
   % RunMADS: Run the MADS optimization algorithm
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     runMADS
   % VARIABLES:
   %  param                        = parameter indicating type of MADS run
   %  GUI.Options.RunUntilFeasible = flag to run MADS only until feasible
   %***********************************************************************
case 'RunMADS'
   if (param == -1), GUI.Options.RunUntilFeasible = 1; end
   runMADS(param);

   %***********************************************************************
   % ViewResults: View the summary results of a MADS run
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % Calls:     viewResults
   % VARIABLES:
   %  GUI               = structure of all GUI handles and variables
   %***********************************************************************
case 'ViewResults'
   viewResults(get(gcbo,'Position'), GUI);

   %***********************************************************************
   % BestX/BestP: View continuous/categorical variables of MADS solutions
   % ----------------------------------------------------------------------
   % Called by: viewResults
   % VARIABLES:
   %  BestLabel = label for viewing window
   %  best      = user data for categorical variables pushbutton
   %  n         = length of best
   %  Add0Pts   = structure of viewing options
   %  Best      = temporary storage of categorical variables
   %***********************************************************************
case 'BestX'
   BestLabel = 'Best_Continuous_Variable_Values';
   assignin('base',BestLabel,get(gcbo,'UserData'));
   openvar(BestLabel);
   clear(BestLabel);
case 'BestP'
   best = get(gcbo,'UserData');   
   BestLabel = 'Best_Categorical_Variable_Values';
   n = length(best);
   for k = 1:n
      if (~ischar(best{k})), best{k} = num2str([best{k}]); end
   end
   Add0pts = struct('Resize','on','WindowStyle','normal', ...
                    'Interpreter','tex');
   Best = inputdlg(cellstr(int2str([1:n]')),BestLabel,ones(n,1)*[1,50], ...
                   best,Add0pts);
   clear('Best');

   %***********************************************************************
   % CopyPlot: Copy history or filter plot to a new screen
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  h1   = handle to new figure window
   %  h2   = handle to new axes on figure window
   %***********************************************************************
case 'CopyPlot'
   h1 = figure;
   h2 = copyobj(gcbo,h1);
   set(h2,'Position','default','ButtonDownFcn','');

   %***********************************************************************
   % Help: View a Help file
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  errorflag = error flag for viewing failure
   %***********************************************************************
case 'Help'
   errorflag = web(['file:///', which(get(gcbo,'UserData'))], '-browser');
   if (errorflag), errordlg('Error: browser or help file not found'); end

   %***********************************************************************
   % About: View author information
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI.Program  =   Software name
   %  GUI.PVersion =   Software version number
   %***********************************************************************
case 'About'
   msgbox({[GUI.Program ',  Version ' GUI.PVersion], ...
          'Copyright (c) 2001-2004 by Mark A. Abramson'}, ...
          ['About ' GUI.Program]);

   %***********************************************************************
   % StopRun: Stop a MADS run
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI.StopRun  =   handle for Stop Run pushbutton
   %***********************************************************************
case 'StopRun'
   set(gcbo,'UserData',1);

   %***********************************************************************
   % Quit: End program session
   % ----------------------------------------------------------------------
   % Called by: NOMADm
   % VARIABLES:
   %  GUI.Path =   path of current optimization problem
   %  GUI.fig  = handle for GUI figure window
   %***********************************************************************
case 'Quit'
   rmpath(GUI.Path);
   if ishandle(GUI.fig)
      delete(GUI.fig);
   end
   clear variables;
   clear functions;
end

drawnow;
return

%**************************************************************************
% UpdateScreen:  Displays a user input screen to set Search parameters.
% -------------------------------------------------------------------------
% Called by: nomadm_functions
% Calls:     updateSearchLabels
% VARIABLES:
%  newGUI              = updated GUI
%  updateType          = type of screen update
%  Choice              = structure of menu choices
%    .PollStrategy     =   integer choice of Poll strategy
%    .PollOrder        =   integer choice of Poll order strategy
%    .PollCenter       =   integer choice of Poll center
%  Options             = structure of MADS parameters
%    .nSearches        =   number of Search types used
%    .Search(n)        =   structure of Search parameters
%    .delta0           =   initial mesh size
%    .deltaMax         =   maximum mesh size
%    .meshRefine       =   mesh refinement factor
%    .meshCoarsen      =   mesh coarsening factor
%    .CacheTol         =   tolerance for flagging point as being in Cache
%    .hmin             =   minimum h-value of an infeasible point
%    .hmax             =   maximum h-value of a filter point
%    .EPollTriggerF    =   f-value Extended Poll trigger
%    .EPollTriggerH    =   h-value Extended Poll trigger
%    .Term             =   substructure containing MADS termination criteria
%      .delta          =     mesh size parameter
%      .iter           =     maximum number of MADS iterations
%      .func           =     maximum number of function evaluations
%      .time           =     maximum CPU time
%      .fails          =     maximum number of consecutive Poll failures
%    .TermIterFlag     =   turns on/off .iter as termination criteria
%    .TermFuncFlag     =   turns on/off .func as termination criteria
%    .TermTimeFlag     =   turns on/off .time as termination criteria
%    .TermFailsFlag    =   turns on/off .fails as termination criteria
%    .scale            =   flag for scaling mesh directions
%    .PollComplete     =   turns on/off complete Polling
%    .TermRel          =   computes termination delta relative to .delta0
%    .useFilter        =   use filter for nonlinear constraints
%    .removeRedundancy = remove redundant linear constraints
%    .accelerate       =   flag for accelerating mesh refinement
%    .plotHistory1     =   turns on/off a history plot
%    .plotHistory2     =   turns on/off a real-time history plot
%    .plotFilter       =   turns on/off a real-time filter plot
%    .loadCache        =   flag for loading a pre-existing Cache file
%    .countCache       =   flag for counting Cache points as function calls
%    .RunOneIteration  =   flag for running one MADS iteration at a time
%    .RunUntilFeasible =   flag for running MADS only until feasible
%  GUI                 = structure for GUI parameters
%                        (These are all handles for GUI window objects)
%    .SearchType(k)    =   current Search types
%    .PollStrategy     =   current Poll strategy
%    .PollCenter       =   current Poll center
%    .PollOrder        =   current Poll order type
%    .delta0           =   current initial mesh size
%    .deltaMax         =   current maximum mesh size
%    .meshRefine       =   current mesh refinement factor
%    .meshCoarsen      =   current mesh coarsening factor
%    .CacheTol         =   current Cache tolerance
%    .hmin             =   current minimum infeasible h-value
%    .hmax             =   current maximum filter h-value
%    .EPollXiF         =   current f-value Extended Poll trigger
%    .EPollXiH         =   current h-value Extended Poll trigger
%    .TermDelta        =   current mesh size termination criteria
%    .TermIter         =   current maximum number of MADS iterations
%    .TermFunc         =   current maximum number of function evaluations
%    .TermTime         =   current maximum CPU time
%    .TermFails        =   current max number of consec Poll failures
%    .TermIterFlag     =   current checkbox value for .TermIter
%    .TermFuncFlag     =   current checkbox value for .TermFunc
%    .TermTimeFlag     =   current checkbox value for .TermTime
%    .TermFailsFlag    =   current checkbox value for .TermFails
%    .CacheMenuXXX     =   Cache menu items
%    .MADSMenuXXX      =   MADS menu items
%    .OptionsMenuXXX   =   Options menu items
%    .RunMenuXXX       =   Run menu items
%  onoff               = cell array of two strings "on" and "off"
%  maxDisplay          = number of Search Types displayed on the main GUI
%  searchLabel         = Search label that appears on the main GUI
%**************************************************************************
function newGUI = updateScreen(updateType,Choice,Options,GUI);

onoff = {'on','off'};

% Update GUI Search fields
maxDisplay = length(GUI.SearchType);
searchLabel = updateSearchLabels(maxDisplay,Options.nSearches,Options.Search);
for k = 1:maxDisplay
   set(GUI.SearchType(k), 'String', searchLabel{k});
end

% Update other GUI fields
if (strcmp(updateType, 'all'))
   set(GUI.PollStrategy, 'String', ...
       GUI.Labels.PollStrategy(Choice.pollStrategy));
   set(GUI.PollCenter,   'String', ...
       GUI.Labels.PollCenter(Choice.pollCenter));
   set(GUI.PollOrder,    'String', ...
       GUI.Labels.PollOrder(Choice.pollOrder));
   set(GUI.delta0,       'String', num2str(Options.delta0,       '%1.1f'));
   set(GUI.deltaMax,     'String', num2str(Options.deltaMax,     '%2.1f'));
   set(GUI.MeshRefine,   'String', num2str(Options.meshRefine,   '%1.1f'));
   set(GUI.MeshCoarsen,  'String', num2str(Options.meshCoarsen,  '%1.1f'));
   set(GUI.CacheTol,     'String', num2str(Options.CacheTol,     '%1.6g'));
   set(GUI.hmin,         'String', num2str(Options.hmin,         '%1.5g'));
   set(GUI.hmax,         'String', num2str(Options.hmax,         '%1.2f'));
   set(GUI.EPollXiF,     'String', num2str(Options.EPollTriggerF,'%1.5g'));
   set(GUI.EPollXiH,     'String', num2str(Options.EPollTriggerH,'%1.5g'));
   set(GUI.TermDelta,    'String', num2str(Options.Term.delta,   '%1.5g'));
   set(GUI.TermIter,     'String', Options.Term.iter);
   set(GUI.TermFunc,     'String', Options.Term.func);
   set(GUI.TermTime,     'String', Options.Term.time);
   set(GUI.TermFails,    'String', Options.Term.fails);
   set(GUI.TermIterFlag,  'Value', Options.TermIterFlag);
   set(GUI.TermFuncFlag,  'Value', Options.TermFuncFlag);
   set(GUI.TermTimeFlag,  'Value', Options.TermTimeFlag);
   set(GUI.TermFailsFlag, 'Value', Options.TermFailsFlag);
   set(GUI.ScaleMenu(1),   'Checked', onoff{1+~(Options.scale == 0)});
   set(GUI.ScaleMenu(2),   'Checked', onoff{1+~(Options.scale == 2)});
   set(GUI.ScaleMenu(3),   'Checked', onoff{1+~(Options.scale == 10)});
   set(GUI.CacheMenuLoad,  'Checked', onoff{2-Options.loadCache});
   set(GUI.CacheMenuCount, 'Checked', onoff{2-Options.countCache});
   set(GUI.MADSMenuPollComplete, ...
                           'Checked', onoff{2-Options.PollComplete});
   set(GUI.OptionsMenuTermRel, ...
                           'Checked', onoff{2-Options.TermRel});
   set(GUI.OptionsMenuUseFilter, ...
                           'Checked', onoff{2-Options.useFilter});
   set(GUI.OptionsMenuRemoveRedundancy, ...
                           'Checked', onoff{2-Options.removeRedundancy});
   set(GUI.OptionsMenuAccelerate, ...
                           'Checked', onoff{2-Options.accelerate});
   set(GUI.OptionsMenuPlotHistory1, ...
                           'Checked', onoff{2-Options.plotHistory1});
   set(GUI.OptionsMenuPlotHistory2, ...
                           'Checked', onoff{2-Options.plotHistory2});
   set(GUI.OptionsMenuPlotFilter, ...
                           'Checked', onoff{2-Options.plotFilter});
   set(GUI.RunMenuOneIteration, ...
                           'Checked', onoff{2-Options.RunOneIteration});
   set(GUI.RunMenuExecFeasible, ...
                           'Checked', onoff{2-Options.RunUntilFeasible});
end
newGUI = GUI;
return

%**************************************************************************
% updateSearchLabels:  Displays the Search Labels on the main GUI.
% -------------------------------------------------------------------------
% Called by: updateScreen
% VARIABLES:
%  label          = cell array containing labels for each Search type
%  maxDisplay     = number of Search types shown on the NOMADm GUI
%  nSearches      = number of Search types
%  Search(k)      = structure that describes the k-th Search
%    .type        =   string identifying Search type
%    .label       =   long text label for Search type
%    .nIter       =   number of iterations to perform Search
%    .nPoints     =   number of Search points
%    .file        =   optional user file defining Search
%  k              = Search type counter
%  searchFilePath = path for optional Search file
%  searchFile     = name of optional Search file
%**************************************************************************
function label = updateSearchLabels(maxDisplay, nSearches, Search)

% Initialize display parameters
[label{1:maxDisplay}] = deal('None');
if (nSearches > maxDisplay)
   label{maxDisplay} = 'Multiple Search Types';
   nSearches = maxDisplay - 1;
end

% Construct Search labels
for k = 1:nSearches
    [searchFilePath,searchFile] = fileparts(Search(k).file);
    switch Search(k).type
    case {'None'}
       label{k} = Search(k).label; 
    case {'LHS', 'Mesh','GA'}
       label{k} = [int2str(Search(k).nPoints),'-pt. ',Search(k).label];
    case {'SPollI','CPollI','GPollI'}
       if (Search(k).nPoints == 1)
         label{k} = strrep(Search(k).label(1:end-1), 'N ', '');
       else
         label{k} = strrep(Search(k).label,'N',int2str(Search(k).nPoints));
       end
    case {'DACE'}
       label{k} = ['DACE Surrogate: ',   searchFile];
    case {'Custom'}
       label{k} = ['Custom Search: ',    searchFile];
    case {'CustomS'}
       label{k} = ['Custom Surrogate: ', searchFile];
    end
    label{k} = [label{k}, ' (', int2str(Search(k).nIter), ')'];
end
return

%**************************************************************************
% runMADS:  Assign GUI input fields to MADS variables and run MADS
% -------------------------------------------------------------------------
% Called by: nomadm_functions
% Calls:     loadMADS, mads, nomad_functions('SaveCache')
% VARIABLES:
%  restart           = flag for resuming previous run
%  GUI               = structure containing all GUI handles and variables
%    .Program        =   Name of this software
%    .RunCount       =   current MADS run number
%    .RunMax         =   maximum allowed MADS run number
%    .fig            =   handle for the main figure window
%    .StopRun        =   handle for Stop Run pushbutton
%    .ResumeRun      =   handle for Resume Run pushbutton
%    .axesFilter     =   handle for Filter plot
%    .axesHistory    =   handle for History plot
%    .RunMenuXXX     =   handles for Run menu items
%    .ResultsMenuXXX =   handles for Results menu items
%    .CacheMenuXXX   =   handles for Cache menu items
%    .RunStatus      =   handle for figure window status bar
%  RunSet            = global variable used to store MADS output data
%    .BestF          =   best feasible solution found by MADS
%    .BestI          =   best infeasible solution found by MADS
%    .RunStats       =   structure of MADS run statistics
%      .delta        =     current mesh size
%      .time         =     current CPU time expended
%      .Cache        =   structure of all processed MADS iterates
%  Problem           = structure containing optimization problem data
%    .File.I         =   name of initial points file
%  Options           = structure containing MADS parameters
%    .loadCache      =   flag for loading a pre-existing Cache file
%    .countCache     =   flag for counting Cache points as function calls
%    .delta0         =   initial mesh size
%    .RunCount       =   current MADS run number
%  iterate0          = initial iterate
%  BestF             = best feasible iterate found by MADS
%  BestI             = least infeasible iterate found by MADS
%  RunStats          = structure containing MADS run statistics
%**************************************************************************
function runMADS(restart);

global GUI
if isappdata(0,'RUNSET'), RunSet = getappdata(0,'RUNSET'); end

lasterr([GUI.Program ' interrupted by user']);
try
   % Change figure window, as appropriate
   if (GUI.RunCount >= GUI.RunMax)
      error('Too many MADS runs without clearing.');
   end
   setptr(GUI.fig,'watch');
   set(GUI.StopRun,  'Visible','on','UserData',0);
   set(GUI.ResumeRun,'Visible','off');
   set(GUI.RunMenuClear,'Enable','on');
   set([GUI.axesFilter; get(GUI.axesFilter,'Children')], 'Visible', 'off');
   set(GUI.axesFilter,'ButtonDownFcn','');
   drawnow;

   % Load MADS input data
   [Problem,Options] = loadMADS(restart, GUI);

   % Get initial point
   if (restart > 0)
      iterate0 = [RunSet(GUI.RunCount).BestF, RunSet(GUI.RunCount).BestI];
   elseif (~restart & exist(Problem.File.I))
      iterate0 = feval(Problem.File.I);
   else
      error(['Cannot find: ', Problem.File.I, '.m']);
   end

   % Set up and run MADS algorithm and store output
   if (restart == 2)
      nomadm_functions('SaveCache');
      set(GUI.axesHistory, 'NextPlot','replacechildren');
      Options.loadCache   = 1;
      Options.countCache = 1;
      Options.delta0 = RunSet(GUI.RunCount).RunStats.delta;
      Options.plotColor = GUI.Types.PlotColors{GUI.RunCount};
      set(GUI.RunStatus,'String', ...
                       ['Resuming Run # ', int2str(GUI.RunCount)]);
      RunStats = RunSet(GUI.RunCount).RunStats;
      [BestF,BestI,RunStats,RunSet(1).Cache] = mads(Problem,iterate0, ...
                                                    Options,RunStats);
      delete(fullfile(GUI.Path,Problem.File.C));
      set([GUI.CacheMenuLoad; GUI.CacheMenuCount; GUI.CacheMenuDelete], ...
          'Enable','off');
   else
      set(GUI.RunStatus,'String', ...
                       ['Run # ', int2str(GUI.RunCount+1), ' processing']);
      [BestF,BestI,RunStats,RunSet(1).Cache] = mads(Problem,iterate0, ...
                                                    Options);
      GUI.RunCount = GUI.RunCount+1;
   end
   RunSet(GUI.RunCount).BestF    = BestF;
   RunSet(GUI.RunCount).BestI    = BestI;
   RunSet(GUI.RunCount).RunStats = RunStats;
   if (RunStats.time > 60)
      load train;
      sound(y);
   end

% Perform these tasks if error in MADS run
catch
   set(GUI.RunStatus,'String', ...
                    ['Run # ', int2str(GUI.RunCount+1), ' failed']);
   set([GUI.axesFilter; get(GUI.axesFilter, 'Children')],'Visible','off');
   if (GUI.RunCount == 0)
      set([GUI.axesHistory; get(GUI.axesHistory,'Children')],'Visible', ...
                                                             'off');
   else
      set(GUI.axesHistory,'ButtonDownFcn', ...
                          'nomadm_functions(''CopyPlot'');');
   end
   RunSet(1).Cache = [];
   errordlg(lasterr,'MADS Runtime Error','modal'); beep
   error(lasterr);
end

% Change figure window, as appropriate
set(GUI.axesHistory, 'NextPlot','add');
set(GUI.axesHistory, 'ButtonDownFcn','nomadm_functions(''CopyPlot'');');
set(GUI.axesFilter,  'ButtonDownFcn','nomadm_functions(''CopyPlot'');');

if (GUI.RunCount > 0)
   set([GUI.CacheMenuSave;GUI.RunMenuResume;GUI.RunMenuRestart],...
       'Enable','on');
   set([GUI.ResultsMenu; GUI.ResultsMenuView(GUI.RunCount)], ...
       'Visible','on');
   set(GUI.RunStatus,'String', ...
                    ['Run # ', int2str(GUI.RunCount), ' complete']);
end
set(GUI.StopRun,   'Visible','off','UserData',0);
set(GUI.ResumeRun, 'Visible','on');
setptr(GUI.fig,'arrow');
setappdata(0,'RUNSET',RunSet);

return

%**************************************************************************
% loadMADS:  Assign GUI input fields to MADS variables and run MADS
% -------------------------------------------------------------------------
% Called by: nomadm_functions('SaveSession'), runMADS
% Calls:     nomadm_compile
% VARIABLES:
%  Problem             = structure containing optimization problem data
%    .File             =   structure of problem file names
%      .F              =   name of functions file
%      .O              =   name of Omega file
%      .I              =   name of initial points file
%      .N              =   name of discrete neighbors file
%      .C              =   name of Cache File
%      .CacheName      =   name of the base workspace Cache variable
%      .FType          =   type of functions file (M=MATLAB,F=FORTRAN,C=C)
%  Options             = structure containing MADS parameters
%    .dacePath         =   path of the DACE Toolbox software
%    .nSearches        =   number of Search types used
%    .Search(n)        =   structure of Search parameters
%      .type           =     string identifying the type of Search
%    .PollStrategy     =   string identifying selected Poll strategy
%    .PollOrder        =   string identifying selected Poll order strategy
%    .PollCenter       =   integer identifying selected Poll center
%    .PollComplete     =   turns on/off complete Polling
%    .loadCache        =   flag for loading a pre-existing Cache file
%    .countCache       =   flag for counting Cache points as function calls
%    .useFilter        =   use filter for nonlinear constraints
%    .removeRedundancy =   remove redundant linear constraints
%    .accelerate       =   flag for accelerating mesh refinement
%    .scale            =   flag for scaling mesh directions
%    .plotHistory1     =   turns on/off a history plot
%    .plotHistory2     =   turns on/off a real-time history plot
%    .plotFilter       =   turns on/off a real-time filter plot
%    .RunOneIteration  =   flag for running one MADS iteration at a time
%    .RunUntilFeasible =   flag for running MADS only until feasible
%    .TermRel          =   computes termination delta relative to .delta0
%    .RunCount         =   MADS run counter
%    .hplothandle      =   handle for history plot axes
%    .fplothandle      =   handle for filter plot axes
%    .stophandle       =   handle for Stop Run pushbutton
%    .delta0           =   initial mesh size
%    .deltaMax         =   maximum mesh size
%    .meshRefine       =   mesh refinement factor
%    .meshCoarsen      =   mesh coarsening factor
%    .CacheTol         =   tolerance for flagging point as being in Cache
%    .hmin             =   minimum h-value of an infeasible point
%    .hmax             =   maximum h-value of a filter point
%    .EPollTriggerF    =   f-value Extended Poll trigger
%    .EPollTriggerH    =   h-value Extended Poll trigger
%    .Term             =   substructure containing MADS termination criteria
%      .delta          =     mesh size parameter
%      .iter           =     maximum number of MADS iterations
%      .func           =     maximum number of function evaluations
%      .time           =     maximum CPU time
%      .fails          =     maximum number of consecutive Poll failures
%    .TermIterFlag     =   turns on/off .iter as termination criteria
%    .TermFuncFlag     =   turns on/off .func as termination criteria
%    .TermTimeFlag     =   turns on/off .time as termination criteria
%    .TermFailsFlag    =   turns on/off .fails as termination criteria
%  restart             = flag for resuming previous run
%  GUI                 = structure containing all GUI handles and variables
%    .Path             =   path of current optimization problem
%    .ProblemExt       =   filename extension of optimization problem 
%    .ProblemName      =   name of current optimization problem
%    .RunStatus        =   handle for GUI figure window status bar
%    .CacheMenuXXX     =   handles for Cache menu items
%    .MADSMenuXXX      =   handles for MADS menu items
%    .OptionsMenuXXX   =   handles for Options menu items
%    .RunMenuXXX       =   handles for Run menu items
%    .Options          =   current Options values
%    .Types            = lists of possible type
%      .Search         =   list of possible Search types
%      .Poll           =   list of possible Poll strategies
%      .PollOrder      =   list of possible Poll order strategies
%    .Choice           = user choices
%      .pollStrategy   =   selected Poll strategy
%      .pollCenter     =   selected Poll center
%      .pollOrder      =   selected Poll order strategy
%                          (These are all handles for GUI window objects)
%    .delta0           =   current initial mesh size
%    .deltaMax         =   current maximum mesh size
%    .meshRefine       =   current mesh refinement factor
%    .meshCoarsen      =   current mesh coarsening factor
%    .CacheTol         =   current Cache tolerance
%    .hmin             =   current minimum infeasible h-value
%    .hmax             =   current maximum filter h-value
%    .EPollXiF         =   current f-value Extended Poll trigger
%    .EPollXiH         =   current h-value Extended Poll trigger
%    .TermDelta        =   current mesh size termination criteria
%    .TermIter         =   current maximum number of MADS iterations
%    .TermFunc         =   current maximum number of function evaluations
%    .TermTime         =   current maximum CPU time
%    .TermFails        =   current max number of consec Poll failures
%    .TermIterFlag     =   current checkbox value for .TermIter
%    .TermFuncFlag     =   current checkbox value for .TermFunc
%    .TermTimeFlag     =   current checkbox value for .TermTime
%    .TermFailsFlag    =   current checkbox value for .TermFails
%  ProblemName         = name of the optimization problem to be solved
%  language            = programming language of functions file
%  k                   = Search counter
%**************************************************************************
function [Problem,Options] = loadMADS(restart,GUI);

% Transfer Optimization Problem data from GUI into MADS input variables
addpath(GUI.Path);
ProblemName    = get(GUI.ProblemName, 'String');
%Problem.File  = [ProblemName, GUI.FileExt];
Problem.File.F = [ProblemName, GUI.FileExt.F];
Problem.File.O = [ProblemName, GUI.FileExt.O];
Problem.File.I = [ProblemName, GUI.FileExt.I];
Problem.File.N = [ProblemName, GUI.FileExt.N];
Problem.File.P = [ProblemName, GUI.FileExt.P];
Problem.File.C = [ProblemName, GUI.FileExt.C];
Problem.CacheName = GUI.CacheName;

% Compile non-Matlab Functions file, if possible
Problem.FType  = upper(GUI.ProblemExt(2));
if ~strcmp(Problem.FType,'M')
   language = nomadm_compile(Problem.FType, ...
                             GUI.Path,Problem.File.F,GUI.ProblemExt);
   set(GUI.RunStatus,'String',['Compiling ',language,' function file']);
end

% Transfer user options from GUI into MADS input variables
umtoggle(GUI.ScaleMenu(2));
umtoggle(GUI.ScaleMenu(3));
umtoggle(GUI.CacheMenuLoad);
umtoggle(GUI.CacheMenuCount);
umtoggle(GUI.MADSMenuPollComplete);
umtoggle(GUI.OptionsMenuUseFilter);
umtoggle(GUI.OptionsMenuRemoveRedundancy);
umtoggle(GUI.OptionsMenuAccelerate);
umtoggle(GUI.OptionsMenuPlotFilter);
umtoggle(GUI.OptionsMenuPlotHistory1);
umtoggle(GUI.OptionsMenuPlotHistory2);
umtoggle(GUI.RunMenuExecFeasible);
umtoggle(GUI.RunMenuOneIteration);
umtoggle(GUI.OptionsMenuTermRel);

Options = GUI.Options;
Options.dacePath = GUI.dacePath;
Options.Search(Options.nSearches+1:end) = [];
for k = 1:Options.nSearches
   Options.Search(k).type = GUI.Types.Search{GUI.Choice.search(k)};
end
Options.PollStrategy     = GUI.Types.Poll{GUI.Choice.pollStrategy};
Options.PollOrder        = GUI.Types.PollOrder{GUI.Choice.pollOrder};
Options.PollCenter       = GUI.Choice.pollCenter - 1;
Options.PollComplete     = umtoggle(GUI.MADSMenuPollComplete);
Options.loadCache        = umtoggle(GUI.CacheMenuLoad);
Options.countCache       = umtoggle(GUI.CacheMenuCount);
Options.useFilter        = umtoggle(GUI.OptionsMenuUseFilter);
Options.removeRedundancy = umtoggle(GUI.OptionsMenuRemoveRedundancy);
Options.accelerate       = umtoggle(GUI.OptionsMenuAccelerate);
Options.scale            =  2*umtoggle(GUI.ScaleMenu(2)) + ...
                           10*umtoggle(GUI.ScaleMenu(3));
Options.plotFilter       = umtoggle(GUI.OptionsMenuPlotFilter);
Options.plotHistory1     = umtoggle(GUI.OptionsMenuPlotHistory1);
Options.plotHistory2     = umtoggle(GUI.OptionsMenuPlotHistory2);
Options.RunOneIteration  = umtoggle(GUI.RunMenuOneIteration);
Options.RunUntilFeasible = umtoggle(GUI.RunMenuExecFeasible);
Options.TermRel          = umtoggle(GUI.OptionsMenuTermRel);
Options.plotColor        = GUI.Types.PlotColors{GUI.RunCount+1};
Options.hplothandle      = GUI.axesHistory;
Options.fplothandle      = GUI.axesFilter;
Options.stophandle       = GUI.StopRun;
Options.delta0           = str2double(get(GUI.delta0,      'String'));
Options.deltaMax         = str2double(get(GUI.deltaMax,    'String'));
Options.meshRefine       = str2double(get(GUI.MeshRefine,  'String'));
Options.meshCoarsen      = str2double(get(GUI.MeshCoarsen, 'String'));
Options.CacheTol         = str2double(get(GUI.CacheTol,    'String'));
Options.hmin             = str2double(get(GUI.hmin,        'String'));
Options.hmax             = str2double(get(GUI.hmax,        'String'));
Options.EPollTriggerF    = str2double(get(GUI.EPollXiF,    'String'));
Options.EPollTriggerH    = str2double(get(GUI.EPollXiH,    'String'));
Options.Term.delta       = str2double(get(GUI.TermDelta,   'String'));
Options.Term.iter        = str2double(get(GUI.TermIter,    'String'));
Options.Term.func        = str2double(get(GUI.TermFunc,    'String'));
Options.Term.time        = str2double(get(GUI.TermTime,    'String'));
Options.Term.fails       = str2double(get(GUI.TermFails,   'String'));
Options.TermIterFlag     = get(GUI.TermIterFlag,  'Value');
Options.TermFuncFlag     = get(GUI.TermFuncFlag,  'Value');
Options.TermTimeFlag     = get(GUI.TermTimeFlag,  'Value');
Options.TermFailsFlag    = get(GUI.TermFailsFlag, 'Value');
if (~Options.TermIterFlag),  Options.Term.iter  = Inf; end
if (~Options.TermFuncFlag),  Options.Term.func  = Inf; end
if (~Options.TermTimeFlag),  Options.Term.time  = Inf; end
if (~Options.TermFailsFlag), Options.Term.fails = Inf; end
if (Options.TermRel)
   Options.Term.delta = Options.Term.delta*Options.delta0;
end
return

%**************************************************************************
% searchOptions:  Displays a user input screen to set Search parameters.
% -------------------------------------------------------------------------
% Called by: nomadm_functions
% VARIABLES:
%  n                 = number of Search types to be used during MADS run
%  optimizer         = string identifying the surrogate optimizer
%  GUI               = structure containing all GUI handles and variables
%    .figSearch      = handle for Search figure window
%    .maxSearches    = maximum number of Search types that can be selected
%    .SearchScreen   = structure of handles for each k of n Searches
%      .type         =   handle for popupmenu of possible Search types
%      .nIter        =   handle for number of iterations field
%      .nPoints      =   handle for number of Search points field
%      .file         =   handle for string containing name of Search file
%    .Labels.Search  = long text labels used in Search type popup menus
%    .Choice.Search  = integers recording Search type popup menu choices
%    .Options.Search = vector of structures for each k of n Searches
%      .nIter        =   number of iterations
%      .nPoints      =   number of Search points
%      .file         =   string containing name of optional Search file
%    .SubsFile       =   name of file that controls GUI callback functions
%  k                 = Search type counter
%  row               = row location of current Search figure window object
%  visible           = flag for making unused Searches invisible
%  PathStr           = temporary storage for Search file path
%  FileName          = temporary storage for Search filename
%  FileExt           = temporary storage for Search filename extension
%**************************************************************************
function searchOptions(n, optimizer);

global GUI

%Set up "Search Options" figure window
GUI.figSearch = figure(...
   'Name', 'Set Options for MADS SEARCH step', ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.375, ...,
   'DefaultUIControlBackgroundColor', [.8 .8 .8], ...
   'DefaultUIControlStyle',           'text', ...
   'Units',                           'normalized', ...
   'Position',                        [0.15 0.1 0.8 0.7], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Figure window shading to make it appear recessed
uicontrol(GUI.figSearch, 'Style','frame','Position', [0, 0, 1, .004],   ...
    'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figSearch, 'Style','frame','Position', [.997, 0 .003, 1], ...
    'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figSearch, 'Style','frame','Position', [0, .996, 1, .004],...
    'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(GUI.figSearch, 'Style','frame','Position', [0, 0, .003, 1],   ...
    'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Additional figure window partitioning lines
uicontrol(GUI.figSearch, 'Style','frame','Position', [0, .896, 1, .004],...
    'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(GUI.figSearch, 'Style','frame','Position', [0, .9,   1, .004],...
    'ForegroundColor','white','BackgroundColor','white');

% Display Number of Search Types and Surrogate Optimizer
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', 'Number of Search Types: ', ...
    'Position', [.13 .905 .20 .06]);
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', int2str(n), ...
    'Position', [.35 .905 .05 .06]);
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', 'Surrogate Optimizer: ', ...
    'Position', [.53 .905 .15 .06]);
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', optimizer, ...
    'Position', [.75 .905 .10 .06]);

% Text headers for the Search parameters
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', 'Search Strategy', ...
    'Position', [.12 .80 .31 .06]);
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', 'Number of Iterations', ...
    'Position', [.45 .80 .15 .06]);
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', 'Number of Points', ...
    'Position', [.62 .80 .16 .06]);
uicontrol(GUI.figSearch, 'Style', 'text', ...
    'String', 'User File', ...
    'Position', [.80 .80 .16 .06]);

% Main loop for each of n Searches
for k = 1:GUI.maxSearches
   if (k <= n)
      row = .76 - .08*(k-1);
      visible = 'on';
   else
      visible = 'off';
   end

   uicontrol(GUI.figSearch, ...
      'Style',           'text', ...
      'String',          ['Search #', int2str(k), ':'], ...
      'Visible',         visible, ...      
      'Position',        [.03, row-.0075, .08, .06]);

   % The data fields for each Search
   GUI.SearchScreen(k).type = uicontrol(GUI.figSearch, ...
      'Style',           'popupmenu', ...
      'String',          GUI.Labels.Search, ...
      'Visible',         visible, ...
      'BackgroundColor', 'white', ...
      'Value',           GUI.Choice.search(k), ...
      'Position',        [.12, row, .31, .06], ...
      'UserData',        k, ...
      'Callback',        [GUI.SubsFile '(''LoadUserSearchFile'');']);
   GUI.SearchScreen(k).nIter = uicontrol(GUI.figSearch, ...
      'Style',           'edit', ...
      'String',          int2str(GUI.Options.Search(k).nIter), ...
      'Visible',         visible, ...
      'BackgroundColor', 'white', ...
      'FontSize',        .64, ...
      'Position',        [.45, row+.02, .15, .04]);
   GUI.SearchScreen(k).nPoints = uicontrol(GUI.figSearch, ...
      'Style',           'edit', ...
      'String',          int2str(GUI.Options.Search(k).nPoints), ...
      'Visible',         visible, ...
      'BackgroundColor', 'white', ...
      'FontSize',        .64, ...
      'Position',        [.62, row+.02, .16, .04]);
   [PathStr,FileName,FileExt] = fileparts(GUI.Options.Search(k).file);
   GUI.SearchScreen(k).file = uicontrol(GUI.figSearch, ...
      'Style',           'text', ...
      'String',          [FileName,FileExt], ...
      'Visible',         visible, ...
      'ForegroundColor', 'red', ...
      'Position',        [.80, row, .16, .06], ...
      'UserData',        PathStr);

   % Set defaults for unseen Searches
  if (k > n)
     set(GUI.SearchScreen(k).type,    'String','None');
     set(GUI.SearchScreen(k).nIter,   'String','1');
     set(GUI.SearchScreen(k).nPoints, 'String','1');
     set(GUI.SearchScreen(k).file,    'String','');
  end
end

% The Done and Cancel Buttons
uicontrol(GUI.figSearch, ...
   'Style',      'pushbutton', ...
   'String',     'Done', ...
   'FontWeight', 'bold', ...
   'Position',   [.33, row-.16, .10, .08], ...
   'UserData',   1, ...
   'Callback',   [GUI.SubsFile '(''LoadSearchOptions'');']);
uicontrol(GUI.figSearch, ...
   'Style',      'pushbutton', ...
   'String',     'Cancel', ...
   'FontWeight', 'bold', ...
   'Position',   [.50, row-.16, .10, .08], ...
   'UserData',   0, ...
   'Callback',   [GUI.SubsFile '(''LoadSearchOptions'');']);

return

%**************************************************************************
% daceOptions:  Displays a user input screen to set DACE Toolbox parameters.
% -------------------------------------------------------------------------
% Called by: nomadm_functions
% VARIABLES:
%  k                         = Search number for this DACE screen
%  GUI                       = global variable holding all GUI variables
%    .SubsFile               =   name of GUI callback functions file
%    .figDACE                =   handle for DACE figure window
%    .Labels                 =   labels for DACE function popup menus
%      .daceRegression       =     labels for DACE regression functions
%      .daceCorrelation      =     labels for DACE correlation functions
%    .Choice                 =   integer choices for DACE functions
%      .daceReg              =     choice for DACE regression function
%      .daceCorr             =     choice for DACE correlation function
%    .daceOptionsRegression  =   popup menu of regression functions
%    .daceOptionsCorrelation =   popup menu of correlation functions
%    .daceOptionsTheta       =   field for estimating correlation parameter
%    .daceOptionsLower       =   filed for entering lower bound for theta
%    .daceOptionsUpper       =   field for entering upper bound for theta
%    .daceOptionsIsotropic   =   checkbox for isotropic correlations
%  Options                   = structure of MADS options
%    .dace                   =   user-chosen DACE Toolbox parameters
%      .theta                =   estimate for theta
%      .lower                =   lower bound for theta
%      .upper                =   upper bound for theta
%      .isotropic            =   flag for isotropic theta
%**************************************************************************
function daceOptions(k);

global GUI

% Set up "DACE Options" figure window
GUI.figDACE = figure(...
   'Name', ['DACE Toolbox Options (Search #', int2str(k),')'], ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.5, ...,
   'DefaultUIControlBackgroundColor', [.8 .8 .8], ...
   'DefaultUIControlStyle',           'text', ...
   'WindowStyle',                     'modal', ...
   'Resize',                          'on', ...
   'Units',                           'normalized', ...
   'Position',                        [0.25 0.35 0.40 0.35], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Figure window shading to make it appear recessed
uicontrol(GUI.figDACE, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figDACE, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figDACE, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(GUI.figDACE, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for DACE screen objects and data fields
uicontrol(GUI.figDACE, ...
      'Style',               'text', ...
      'String',              'Regression Model: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .855 .34 .09]);
GUI.daceOptionsReg         = uicontrol(GUI.figDACE, ...
      'Style',               'popupmenu', ...
      'String',              GUI.Labels.daceRegression, ...
      'BackgroundColor',     'white', ...
      'Value',               GUI.Choice.daceReg(k), ...
      'Position',            [.37 .85 .60 .11]);

uicontrol(GUI.figDACE, ...
      'Style',               'text', ...
      'String',              'Correlation Model: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .725 .34 .09]);
GUI.daceOptionsCorr        = uicontrol(GUI.figDACE, ...
      'Style',               'popupmenu', ...
      'String',              GUI.Labels.daceCorrelation, ...
      'BackgroundColor',     'white', ...
      'Value',               GUI.Choice.daceCorr(k), ...
      'Position',            [.37 .72 .60 .11]);

uicontrol(GUI.figDACE, ...
      'Style',               'text', ...
      'String',              'Theta: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .595 .34 .09]);
GUI.daceOptionsTheta       = uicontrol(GUI.figDACE, ...
      'Style',               'edit', ...
      'String',              num2str(GUI.Options.dace(k).Theta), ...
      'FontSize',            .6, ...
      'BackgroundColor',     'white', ...
      'Position',            [.37 .61 .35 .09]);

uicontrol(GUI.figDACE, ...
      'Style',               'text', ...
      'String',              'Lower: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .465 .34 .09]);
GUI.daceOptionsLower       = uicontrol(GUI.figDACE, ...
      'Style',               'edit', ...
      'String',              num2str(GUI.Options.dace(k).Lower), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .48 .35 .09]);

uicontrol(GUI.figDACE, ...
      'Style',               'text', ...
      'String',              'Upper: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .335 .34 .09]);
GUI.daceOptionsUpper       = uicontrol(GUI.figDACE, ...
      'Style',               'edit', ...
      'String',              num2str(GUI.Options.dace(k).Upper), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .35 .35 .09]);

GUI.daceOptionsIsotropic   = uicontrol(GUI.figDACE, ...
      'Style',               'checkbox',  ...
      'String',              'Isotropic', ...
      'Value',               GUI.Options.dace(k).Isotropic,   ...
      'Position',            [.37 .23 .35 .11]);
  
% The Done and Cancel Buttons
uicontrol(GUI.figDACE, ...
   'Style',      'pushbutton', ...
   'String',     'Done', ...
   'FontWeight', 'bold', ...
   'Position',   [.27, .05, .20, .14], ...
   'UserData',   k, ...
   'Callback',   [GUI.SubsFile '(''LoadDACEOptions'');']);
uicontrol(GUI.figDACE, ...
   'Style',      'pushbutton', ...
   'String',     'Cancel', ...
   'FontWeight', 'bold', ...
   'Position',   [.54, .05, .20, .14], ...
   'UserData',   0, ...
   'Callback',   [GUI.SubsFile '(''LoadDACEOptions'');']);

return

%**************************************************************************
% viewResults:  Displays Results from MADS run.
% -------------------------------------------------------------------------
% Called by: nomadm_functions
% VARIABLES:
%  k              = run number
%  GUI            = structure containing handles for the user interface
%    .SubsFile    =   name of GUI callback functions file
%    .figResults  =   parent handle for the MADS Results screen
%    .BestFFValue =   display of BFP objective function value
%    .BestFHValue =   display of BFP constraint violation function value
%    .ViewFXOpt   =   pushbutton for viewing BFP continuous variables
%    .ViewFPOpt   =   pushbutton for viewing BFP categorical variables
%    .BestIFValue =   display of LIP objective function value
%    .BestIHValue =   display of LIP constraint violation function value
%    .ViewIXOpt   =   pushbutton for viewing LIP continuous variables
%    .ViewIPOpt   =   pushbutton for viewing LIP categorical variables
%    .RunDelta    =   display of final mesh size
%    .RunIter     =   display of number of iterations
%    .RunFunc     =   display of number of function evaluations
%    .RunGrad     =   display of number of gradient evaluations
%    .RunTime     =   display of CPU time expended
%    .RunFails    =   display of number of consecutive Poll failures
%    .RunStop     =   display of whether or not run was interrupted by user
%    .RunCHits    =   display of number of Cache hits
%    .TermDelta   =   display of mesh size termination criteria
%    .TermIter    =   display of maximum number of MADS iterations
%    .TermFunc    =   display of maximum number of function evaluations
%    .TermTime    =   display of maximum CPU time
%    .TermFails   =   display of max number of consecutive Poll failures
%  RunSet         = structure containing all information about a MADS run
%    .BestF       =   best feasible iterate for Run k
%    .BestI       =   least infeasible iterate for Run k
%    .RunStats    =   MADS Run k statistics (iterations, CPU time, etc.)
%      .delta     =     final mesh size
%      .iter      =     total number of iterations
%      .func      =     total number of function evaluations
%      .time      =     total CPU time expended
%      .fails     =     final number of consecutive Poll failures
%      .stopRun   =     final value of Stop Run button
%      .CacheHits =     total number of Cache hits
%  digits         = number of digits accuracy to present on screen
%  NoYes          = cell array containing the strings "No" and "Yes"
%  figName        = title of the MADS Results screen
%**************************************************************************
function viewResults(k, GUI);

RunSet = getappdata(0,'RUNSET');
digits = 10;
NoYes = {'No','Yes'};
figName = ['MADS Results for Problem: ', ...
          get(GUI.ProblemName,'String'), ', Run #' int2str(k)];

% Set up "Display Results" figure window
GUI.figResults = figure( ...
   'Name',                           figName, ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.375, ...,
   'DefaultUIControlBackgroundColor', [.8 .8 .8], ...
   'DefaultUIControlStyle',           'text', ...
   'Units',                           'normalized', ...
   'Position',                        [0.1 0.15 0.8 0.7], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Figure window shading to make it appear recessed
uicontrol(GUI.figResults, 'Style','frame','Position', [0, 0, 1, .004], ...
   'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figResults, 'Style','frame','Position', [.997, 0, .003, 1],...
   'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figResults, 'Style','frame','Position', [0, .996, 1, .004],...
   'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(GUI.figResults, 'Style','frame','Position', [0, 0, .003, 1], ...
   'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

uicontrol(GUI.figResults, 'Style','frame','Position', [0, .536, 1, .004],...
   'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(GUI.figResults, 'Style','frame','Position', [0, .54, 1, .004],...
   'ForegroundColor','white','BackgroundColor','white');
uicontrol(GUI.figResults, 'Style','frame','Position', [.5, .54, .003, .46], ...
   'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(GUI.figResults, 'Style','frame','Position', [.497, .54, .003, .46], ...
   'ForegroundColor','white','BackgroundColor','white');

% Set up display of Best Feasible Solution
uicontrol(GUI.figResults, 'String','BEST FEASIBLE SOLUTION', ...
   'FontWeight','bold', ...
   'HorizontalAlignment','center','Position',[.02 .86 .46 .10]);

if (isempty(RunSet(k).BestF))
   uicontrol(GUI.figResults, 'String','NO FEASIBLE ITERATES FOUND', ...
      'HorizontalAlignment','center','Position',[.02 .77 .46 .08]);
else
   uicontrol(GUI.figResults, 'String','Objective Function Value: ', ...
      'HorizontalAlignment','left',  'Position',[.02 .80 .28 .08]);
   uicontrol(GUI.figResults, 'String','Constraint Violation Measure: ', ...
      'HorizontalAlignment','left',  'Position',[.02 .74 .28 .08]);
   GUI.BestFFValue = uicontrol(GUI.figResults, ...
      'String',num2str(RunSet(k).BestF.f,digits), ...
      'HorizontalAlignment','right', 'Position',[.30 .80 .18 .08]);
   GUI.BestFHValue = uicontrol(GUI.figResults, ...
      'String',num2str(RunSet(k).BestF.h,digits), ...
      'HorizontalAlignment','right', 'Position',[.30 .74 .18 .08]);
   GUI.ViewFXOpt = uicontrol(GUI.figResults, ...
      'Style',               'pushbutton', ...
      'String',              'View Continuous Variables', ...
      'FontWeight',          'bold', ...
      'HorizontalAlignment', 'center', ...
      'Position',            [.02 .66 .46 .08], ...
      'UserData',            RunSet(k).BestF.x, ...
      'Callback',            'nomadm_functions(''BestX'');');
   if (~isempty(RunSet(k).BestF.p))
      GUI.ViewFPOpt = uicontrol(GUI.figResults, ...
         'Style',               'pushbutton', ...
         'String',              'View Categorical Variables', ...
         'FontWeight',          'bold', ...
         'HorizontalAlignment', 'center', ...
         'Position',            [.02 .58 .46 .08], ...
         'UserData',            RunSet(k).BestF.p, ...
         'Callback',            'nomadm_functions(''BestP'');');
   end
end

% Set up display of Least Infeasible Solution
uicontrol(GUI.figResults, 'String','LEAST INFEASIBLE SOLUTION', ...
   'FontWeight','bold', ...
   'HorizontalAlignment','center','Position',[.52 .86 .46 .10]);

if (isempty(RunSet(k).BestI))
   uicontrol(GUI.figResults, 'String','NO INFEASIBLE ITERATES FOUND ', ...
      'HorizontalAlignment','center','Position',[.52 .77 .46 .08]);
else
   uicontrol(GUI.figResults, 'String','Objective Function Value: ', ...
      'HorizontalAlignment','left',  'Position',[.52 .80 .28 .08]);
   uicontrol(GUI.figResults, 'String','Constraint Violation Measure: ', ...
      'HorizontalAlignment','left',  'Position',[.52 .74 .28 .08]);
   GUI.BestIFValue = uicontrol(GUI.figResults, ...
      'String',num2str(RunSet(k).BestI.f,digits), ...
      'HorizontalAlignment','right', 'Position',[.80 .80 .18 .08]);
   GUI.BestIHValue = uicontrol(GUI.figResults, ...
      'String',num2str(RunSet(k).BestI.h,digits), ...
      'HorizontalAlignment','right', 'Position',[.80 .74 .18 .08]);

% Set up pushbuttons for viewing solution vectors
   GUI.ViewIXOpt = uicontrol(GUI.figResults, ...
      'Style',               'pushbutton', ...
      'String',              'View Continuous Variables', ...
      'FontWeight',          'bold', ...
      'HorizontalAlignment', 'center', ...
      'Position',            [.52 .66 .46 .08], ...
      'UserData',            RunSet(k).BestI.x, ...
      'Callback',            'nomadm_functions(''BestX'');');
   if (~isempty(RunSet(k).BestI.p))
      GUI.ViewIPOpt = uicontrol(GUI.figResults, ...
         'Style',               'pushbutton', ...
         'String',              'View Categorical Variables', ...
         'FontWeight',          'bold', ...
         'HorizontalAlignment', 'center', ...
         'Position',            [.52 .58 .46 .08], ...
         'UserData',            RunSet(k).BestI.p, ...
         'Callback',            'nomadm_functions(''BestP'');');
   end
end

% Display Labels for MADS Run Statistics
uicontrol(GUI.figResults, 'String','RUN STATISTICS', ...
   'FontWeight','bold', ...
   'HorizontalAlignment','center','Position',[.30 .42 .40 .10]);
uicontrol(GUI.figResults, 'String','Final Mesh Size:', ...
   'HorizontalAlignment','left',  'Position',[.25 .36 .25 .08]);
uicontrol(GUI.figResults, 'String','MADS Iterations:', ...
   'HorizontalAlignment','left',  'Position',[.25 .31 .25 .08]);
uicontrol(GUI.figResults, 'String','Function Evaluations:', ...
   'HorizontalAlignment','left',  'Position',[.25 .26 .25 .08]);
uicontrol(GUI.figResults, 'String','Gradient Evaluations:', ...
   'HorizontalAlignment','left',  'Position',[.25 .21 .25 .08]);
uicontrol(GUI.figResults, 'String','CPU Time:', ...
   'HorizontalAlignment','left',  'Position',[.25 .16 .25 .08]);
uicontrol(GUI.figResults, 'String','Consecutive Poll Failures:', ...
   'HorizontalAlignment','left',  'Position',[.25 .11 .25 .08]);
uicontrol(GUI.figResults, 'String','Interrupted by User:', ...
   'HorizontalAlignment','left',  'Position',[.25 .06 .25 .08]);
uicontrol(GUI.figResults, 'String','Cache Hits:', ...
   'HorizontalAlignment','left',  'Position',[.25 .01 .25 .08]);

% Display MADS Run Statistics (Compare with termination criteria)
GUI.RunDelta = uicontrol(GUI.figResults, ...
   'String',num2str(RunSet(k).RunStats.delta,digits), ...
   'HorizontalAlignment','right',  'Position',[.50 .36 .20 .08]);
GUI.RunIter  = uicontrol(GUI.figResults, ...
   'String',int2str(RunSet(k).RunStats.iter), ...
   'HorizontalAlignment','right',  'Position',[.50 .31 .20 .08]);
GUI.RunFunc  = uicontrol(GUI.figResults, ...
   'String',int2str(RunSet(k).RunStats.func), ...
   'HorizontalAlignment','right',  'Position',[.50 .26 .20 .08]);
GUI.RunGrad  = uicontrol(GUI.figResults, ...
   'String',int2str(RunSet(k).RunStats.grad), ...
   'HorizontalAlignment','right',  'Position',[.50 .21 .20 .08]);
GUI.RunTime  = uicontrol(GUI.figResults, ...
   'String',num2str(RunSet(k).RunStats.time), ...
   'HorizontalAlignment','right',  'Position',[.50 .16 .20 .08]);
GUI.RunFails = uicontrol(GUI.figResults, ...
   'String',num2str(RunSet(k).RunStats.fails), ...
   'HorizontalAlignment','right',  'Position',[.50 .11 .20 .08]);
GUI.RunStop  = uicontrol(GUI.figResults, ...
   'String',NoYes{1+RunSet(k).RunStats.stopRun}, ...
   'HorizontalAlignment','right',  'Position',[.50 .06 .20 .08]);
GUI.RunCHits = uicontrol(GUI.figResults, ...
   'String',num2str(RunSet(k).RunStats.CacheHits), ...
   'HorizontalAlignment','right',  'Position',[.50 .01 .20 .08]);

% Change colors for violated Termination Criteria
if (RunSet(k).RunStats.delta < str2num(get(GUI.TermDelta,'String')))
   set(GUI.RunDelta,'FontWeight','bold','ForegroundColor','blue');
end
if (RunSet(k).RunStats.iter >= str2num(get(GUI.TermIter,'String')))
   set(GUI.RunIter, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.func >= str2num(get(GUI.TermFunc,'String')))
   set(GUI.RunFunc, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.time >= str2num(get(GUI.TermTime,'String')))
   set(GUI.RunTime, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.fails >= str2num(get(GUI.TermFails,'String')))
   set(GUI.RunFails, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.stopRun)
   set(GUI.RunStop, 'FontWeight','bold','ForegroundColor','red');
end
return
