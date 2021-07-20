function nomadm
%NOMADM   Execute the NOMADm graphic user interface (GUI).
%
%   Syntax:
%      nomadm
%
%   Description:
%      NOMADM launches the NOMADm GUI, which controls the setup of an
%      optimization problem, setting of various algorithm parameters and
%      options, running of the MADS algorithm for numerically solving the
%      optimization problem, and viewing results.
%
%   See also MADS, MADS_BATCH, MADS_DEFAULTS, NOMADM_COMPILE

%*******************************************************************************
%   Copyright (c) 2001-2005 by Mark A. Abramson
%
%   This file is part of the NOMADm software package.
%
%   NOMADm is free software; you can redistribute it and/or modify it under the
%   terms of the GNU General Public License as published by the Free Software
%   Foundation; either version 2 of the License, or (at your option) any later
%   version.
%
%   NOMADm is distributed in the hope that it will be useful, but WITHOUT ANY
%   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
%   details.
%
%   You should have received a copy of the GNU General Public License along with
%   NOMADm; if not, write to the Free Software Foundation, Inc., 
%   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% ------------------------------------------------------------------------------
%   Originally created, 2001.
%   Last modified, 2 February 2005
%
%   Author information:
%   Mark A. Abramson, LtCol, USAF, PhD
%   Air Force Institute of Technology
%   Department of Mathematics and Statistics
%   2950 Hobson Way
%   Wright-Patterson AFB, OH 45433
%   (937) 255-3636 x4524
%   Mark.Abramson@afit.edu
%*******************************************************************************

%*******************************************************************************
% nomadm: Runs the NOMADm user interface and associated callback functions.
% ------------------------------------------------------------------------------
% CALLED FUNCTIONS:
%  nomadm_gui             = Launches the NOMADm GUI.
%    selectProblem        =   Choose an optimization problem (CB)
%    editFile             =   Edit a user file (CB)
%    quitGUI              =   Close the GUI and quit the program (CB)
%    editParameters       =   Change certain MADS parameter values (CB)
%    selectPoll           =   Select Poll strategy, center, or order (CB)
%    loadSession          =   Load user options from file (CB)
%    resetSession         =   Reset parameters to default values (CB)
%    saveSession          =   Save user options to file (CB)
%    saveCache            =   Save the Cache to file for later use (CB)
%    deleteFile           =   Delete a Cache or Session file (CB)
%    clearRuns            =   Clears all data from with previous runs (CB)
%    getHelp              =   View a help file (CB)
%    getAbout             =   View information about this software (CB)
%    copyPlot             =   Copy a plot to its own window (BD)
%    runMADS              =   Run the MADS optimizer (CB)
%      loadMADS           =     Transfer GUI values to MADS variables
%    search_gui           =   Launch the Search screen GUI
%      modifySearchScreen =     Change appearance of the Search screen (CB)
%      loadUserSearchFile =     Load a user-specified Search file (CB)
%      loadSearchOptions  =     Load user Search step choices (CB)
%    dace_gui             =   Launch the DACE figure window
%      loadDACEOptions    =     Load user DACE options (CB)
%    nw_gui               =   Launch the NW figure window (CB)
%      loadNWOptions      =     Load user NW options (CB)
%    updateScreen         =   Update contents of the NOMADm figure window
%      updateSearchLabels =     Update Search labels on the NOMADm window
%    results_gui          =   Launch Results GUI, view results of run (CB)
%      solution_gui       =     Launch Solution GUI, view solution (CB)
%        changeView       =       Change view of page with Previous/Next buttons
% ------------------------------------------------------------------------------
% VARIABLES (only for the nomadm function):
%  h                  = handle of pre-existing NOMADm figure window
%  gui_var            = structure of all GUI variables
%    .version         =   Version number for this software
%    .path            =   Current path of the optimization problem files
%    .newPath         =   flag for new path not already on the Matlab path
%    .runCount        =   current MADS run number
%    .runMax          =   maximum allowed MADS run number
%    .noGrad          =   flag indicating initial availability of gradients
%    .Defaults        =   substructure of parameter default values
%    .Labels          =   substructure of NOMADm menu labels
%    .Types           =   substructure of NOMADm types 
%    .FileExt         =   substructure of possible file extensions
%    .HelpDoc         =   substructure of help file names
%    .maxSearches     =   maximum allowable number of Search types used
%    .nameCache       =   name of the Cache that is saved as appdata
%  guiFunction        = cell array of NOMADm function handle names
%  gui_func           = sructure of handles for all NOMADm functions
%  gui                = structure of GUI object handles
%    .fig             =   handle for the main figure window
%*******************************************************************************

% Before launching the NOMADm GUI, make sure it is not already running
h = findobj(0,'type','Figure','Tag','NOMADm_GUI');
if ~isempty(h)
   disp('Only one NOMADm window may run at a time');
   figure(h);
   return
end

% Declare basic GUI variables
gui_var.newMatlab   = str2num(version('-release')) > 12;
gui_var.version     = '3.3';
gui_var.path        = '';
gui_var.newPath     = 0;
gui_var.runCount    = 0;
gui_var.runMax      = 10;
gui_var.noGrad      = 1;

% Set up default values
gui_var.Defaults    = mads_defaults('Truth');
gui_var.Labels      = gui_var.Defaults.Labels;
gui_var.Types       = gui_var.Defaults.Types;
gui_var.FileExt     = gui_var.Defaults.FileExt;
gui_var.HelpDoc     = gui_var.Defaults.HelpDoc;
gui_var.maxSearches = gui_var.Defaults.maxSearches;
gui_var.nameCache   = gui_var.Defaults.nameCache;

% Set function handles for all nomadm functions used in the GUI
guiFunction = {'selectProblem','editFile','quitGUI','editParameters',   ...
               'selectPoll','toggleStochastic','loadSession','resetSession',...
               'saveSession','saveCache', 'deleteFile','clearRuns','getHelp',...
               'getAbout','search_gui','dace_gui','nw_gui','results_gui',...
               'solution_gui','loadSearchOptions','loadDACEOptions',...
               'loadNWOptions','copyPlot','runMADS','loadMADS',...
               'loadUserSearchFile','changeView','updateScreen',...
               'modifySearchScreen','updateSearchLabels'};
gui_func = [];
for k = 1:length(guiFunction)
   gui_func.(guiFunction{k}) = str2func(guiFunction{k});
end

% Launch the NOMADm GUI and store GUI handles and variables as appdata
gui = nomadm_gui(gui_var,gui_func);
setappdata(gui.fig,'gui',gui);
setappdata(gui.fig,'gui_func',gui_func);
setappdata(gui.fig,'gui_var',gui_var);
resetSession(1,0);
return

%*******************************************************************************
% nomadm_gui:  Launch the NOMADm GUI.
% ------------------------------------------------------------------------------
% Calls: mads_defaults
% VARIABLES:
%  gui_var            = structure of GUI variables (descriptions above)
%  gui_func           = structure of function handles with fields as names
%  gui                = structure of all GUI object handles
%    .fig             =   handle for the main figure window
%    .TermFlag        =   structure of handles for termination checkboxes
%      .nIter         =     handle for number of iterations
%      .nFunc         =     handle for number of function evaluations
%      .time          =     handle for CPU time
%      .nFails        =     handle for number of consecutive Poll failures
%    .problem         =   handle for text listing optimization problem name
%    .searchLabel(k)  =   handles for Searches
%    .pollStrategy    =   handle for Poll strategy
%    .pollCenter      =   handle for Poll center
%    .pollOrder       =   handle for Poll order type
%    .Term            =   structure of handles for termination criteria
%      .delta         =     handle for mesh size tolerance
%      .nIter         =     handle for maximum number of MADS iterations
%      .nFunc         =     handle for maximum number of function evaluations
%      .time          =     handle for maximum CPU time
%      .nFails        =     handle for max number of consec Poll failures
%    .delta0          =   handle for initial mesh size
%    .deltaMax        =   handle for maximum mesh size
%    .meshRefine      =   handle for mesh refinement factor
%    .meshCoarsen     =   handle for mesh coarsening factor
%    .tolCache        =   handle for Cache tolerance
%    .hmin            =   handle for minimum infeasible h-value
%    .hmax            =   handle for maximum filter h-value
%    .ePollXiF        =   handle for f-value Extended Poll trigger
%    .ePollXiH        =   handle for h-value Extended Poll trigger
%    .RS              =   structure of handles for ranking & selection
%      .s0            =     handle for initial sample size
%      .alpha_const   =     handle for the initial alpha value
%      .iz_const      =     handle for the initial indifference zone parameter
%      .alpha_rho     =     handle for the alpha decay factor
%      .iz_rho        =     handle for the indifference zone decay factor
%    .runStatus       =   handle for GUI figure window status bar
%    .axesHistory     =   handle for the History plot axes
%    .axesFilter      =   handle for the Filter plot axes
%    .stopRun         =   handle for the Stop Run pushbutton
%    .resumeRun       =   handle for the Resume Run pushbutton
%    .Menu            =   structure of handles for the main menu
%      .Problem       =     handle for the Problem menu
%      .MADS          =     handle for the MADS menu
%      .Options       =     handle for the Options menu
%      .Session       =     handle for the Session menu
%      .Cache         =     handle for the Cache menu
%      .Run           =     handle for the Run menu
%      .Results       =     handle for the Results menu
%      .Help          =     handle for the Help menu
%    .ProblemMenu     =   structure of handles for Problem menu items
%    .MADSMenu        =   structure of handles for MADS menu items
%    .OptionsMenu     =   structure of handles for Options menu items
%    .SessionMenu     =   structure of handles for Session menu items
%    .CacheMenu       =   structure of handles for Cache menu items
%    .RunMenu         =   structure of handles for Run menu items
%    .ResultsMenu     =   structure of handles for Results menu items
%    .HelpMenu        =   structure of handles for Help menu items
%   onoff             = cell array containing "on" and "off" strings
%
%   MANY other variables, primarily ones that are GUI.fig objects
%*******************************************************************************
function gui = nomadm_gui(gui_var,gui_func)
onoff = {'on','off'};

% Set up main figure window
gui.fig = figure(...
   'Name',                               'NOMADm Optimization Software',...
   'Tag',                                'NOMADm_GUI', ...
   'DefaultUIControlUnits',              'normalized', ...
   'DefaultUIControlFontUnits',          'normalized', ...
   'DefaultUIControlFontName',           'Helvetica',  ...
   'DefaultUIControlFontSize',            0.72, ...,
   'DefaultUIControlStyle',              'text', ...
   'DefaultUIControlHorizontalAlignment','left', ...
   'Units',                              'normalized', ...
   'Position',                           [0.05 0.1 .9 .8], ...
   'MenuBar',                            'none', ...
   'NumberTitle',                        'off', ...
   'CloseRequestFcn',                    {gui_func.quitGUI, gui_var});

if gui_var.newMatlab

   % Status bar
   gui.statusbar = uipanel('Parent',gui.fig,'Position',[0, 0, 1, .04], ...
                                            'BorderType','beveledout');
   gui.runStatus = uicontrol(gui.statusbar, 'String','No problem selected', ...
                                            'Position', [.01 0 .98 .88]);

   % Panels for MADS Parameters
   gui.MADSPanel   = uipanel('Parent',gui.fig,'Position',[.02,.745,.395,.20], ...
                                              'BorderType','beveledout');
   gui.TermPanel   = uipanel('Parent',gui.fig,'Position',[.02,.585,.395,.16], ...
                                              'BorderType','beveledout');
   gui.MeshPanel   = uipanel('Parent',gui.fig,'Position',[.02,.425,.395,.16], ...
                                              'BorderType','beveledout');
   gui.FilterPanel = uipanel('Parent',gui.fig,'Position',[.02,.295,.395,.13], ...
                                              'BorderType','beveledout');
   gui.RSPanel     = uipanel('Parent',gui.fig,'Position',[.02,.135,.395,.16], ...
                                              'BorderType','beveledout');
else
   panelColor = [0.9255, 0.9137, 0.8471];
   black = [.35 .35 .35];
   white = [.95 .95 .95];
   
   uicontrol(gui.fig, 'Style','frame','Position',[0, .04, 1, .003], ...
      'ForegroundColor',white,'BackgroundColor',white);
   uicontrol(gui.fig, 'Style','frame','Position',[0,  0,  1, .04], ...
      'ForegroundColor',panelColor,'BackgroundColor',panelColor);
   gui.runStatus = uicontrol(gui.fig, 'String','No problem selected', ...
                                      'Position', [.01 .003 .98 .03], ...
                                      'BackgroundColor',panelColor);

   % Frame box for the Display of MADS Parameters
   gui.MADSPanel   = uicontrol(gui.fig,'Style','frame', ...
                     'Position',       [.02,.745,.395,.20], ...
                     'ForegroundColor',panelColor);
   gui.TermPanel   = uicontrol(gui.fig,'Style','frame', ...
                     'Position',       [.02,.585,.395,.16], ...
                     'ForegroundColor',panelColor);
   gui.MeshPanel   = uicontrol(gui.fig,'Style','frame', ...
                     'Position',       [.02,.425,.395,.16], ...
                     'ForegroundColor',panelColor);
   gui.FilterPanel = uicontrol(gui.fig,'Style','frame', ...
                     'Position',       [.02,.295,.395,.13], ...
                     'ForegroundColor',panelColor);
   gui.RSPanel = uicontrol(gui.fig,'Style','frame',...
                     'Position',       [.02,.135,.395,.16], ...
                     'ForegroundColor',panelColor);
                  
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .295 .004 .65], ...
        'ForegroundColor',white, 'BackgroundColor',white);     % Left side
    uicontrol(gui.fig, 'Style','frame','Position',[.416 .295 .004 .65], ...
        'ForegroundColor',black, 'BackgroundColor',black);     % Right side
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .943 .404 .003], ...
        'ForegroundColor',white, 'BackgroundColor',white);     % Top
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .295 .404 .003], ...
        'ForegroundColor',black, 'BackgroundColor',black);     % Bottom
     
    gui.RSBorder(1) = uicontrol(gui.fig, 'Style','frame',...
        'Position', [.016 .135 .004 .16], ...
        'ForegroundColor',white, 'BackgroundColor',white);     % Left side
    gui.RSBorder(2) = uicontrol(gui.fig, 'Style','frame',...
        'Position', [.416 .135 .004 .16], ...
        'ForegroundColor',black, 'BackgroundColor',black);     % Right side
    gui.RSBorder(3) = uicontrol(gui.fig, 'Style','frame',...
        'Position', [.016 .292 .404 .003], ...
        'ForegroundColor',white, 'BackgroundColor',white);     % Top
    gui.RSBorder(4) = uicontrol(gui.fig, 'Style','frame',...
        'Position', [.016 .135 .404 .003], ...
        'ForegroundColor',black, 'BackgroundColor',black);     % Bottom
 
    % Three separators in the MADS parameter box (2 per separator)
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .745 .404  .003], ...
        'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .742 .404  .003], ...
        'ForegroundColor','white','BackgroundColor','white');
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .585 .404  .003], ...
        'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .582 .404  .003], ...
        'ForegroundColor','white','BackgroundColor','white');
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .425 .404  .003], ...
        'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
    uicontrol(gui.fig, 'Style','frame','Position',[.016 .422 .404  .003], ...
        'ForegroundColor','white','BackgroundColor','white');

end
                                           
% Display of the NOMADm object labels
uicontrol(gui.fig, 'String','Optimization Problem: ', ...
          'FontWeight','bold', ...
          'BackgroundColor',[.8 .8 .8], ...
          'HorizontalAlignment','right', ...
          'Position',[.025 .955  .19  .035], ...
          'ToolTipString',...
          'The name of the currently selected optimization problem');
uicontrol(gui.fig, 'String','MADS Parameter Settings', ...
          'FontWeight','bold', ...
          'HorizontalAlignment','center', ...
          'Position',[.025 .905 .385 .03], ...
          'ToolTipString',...
          'Any of these values may be changed via the MADS Menu');
uicontrol(gui.fig, 'String','Initial Search:', ...
          'Position',[.025 .87  .11  .03], ...
          'ToolTipString','The Strategy used in the initial MADS Search step');
uicontrol(gui.fig, 'String','Search:', ...
          'Position',[.025 .84  .11  .03], ...
          'ToolTipString','The Strategy used in the MADS Search step');
uicontrol(gui.fig, 'String','Poll Directions:', ...
          'Position',[.025 .81  .11  .03], ...
          'ToolTipString',...
          'MADS Poll directions must positively span the problem domain');
uicontrol(gui.fig, 'String','Poll Order:', ...
          'Position',[.025 .78  .11  .03], ...
          'ToolTipString','The order in which the MADS Poll set is evaluated');
uicontrol(gui.fig, 'String','Poll Center:', ...
          'Position',[.025 .75  .11  .03], ...
          'ToolTipString', ...
          'The point around which the MADS Poll step is performed');

gui.TermFlag.delta  = uicontrol(gui.fig, ...
          'String',   gui_var.Labels.Parameters.term{1}, ...
          'Style',    'checkbox', ...
          'Value',    1, ...
          'Position', [.025 .71  .32  .03]);
gui.TermFlag.nIter  = uicontrol(gui.fig, ...
          'String',   gui_var.Labels.Parameters.term{2}, ...
          'Style',    'checkbox', ...
          'Position', [.025 .68 .32 .03]);
gui.TermFlag.nFunc  = uicontrol(gui.fig, ...
          'String',   gui_var.Labels.Parameters.term{3}, ...
          'Style',    'checkbox', ...
          'Position', [.025 .65 .32 .03]);
gui.TermFlag.time   = uicontrol(gui.fig, ...
          'String',   gui_var.Labels.Parameters.term{4}, ...
          'Style',    'checkbox', ...
          'Position', [.025 .62 .32 .03]);
gui.TermFlag.nFails = uicontrol(gui.fig, ...
          'String',   gui_var.Labels.Parameters.term{5}, ...
          'Style',    'checkbox', ...
          'Position', [.025 .59 .32 .03]);

uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{1}, ...
          'Position',[.025 .55  .32  .03], ...
          'ToolTipString','The initial mesh size parameter value');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{2}, ...
          'Position',[.025 .52  .32  .03], ...
          'ToolTipString','The maximum allowed mesh size');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{3}, ...
          'Position',[.025 .49  .32  .03], ...
          'ToolTipString', ...
          'Mesh size is reduced by this factor when an iteration fails');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{4}, ...
          'Position',[.025 .46  .32  .03], ...
          'ToolTipString',...
          'Mesh size is increased by this factor when an iteration succeeds');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.mesh{5}, ...
          'Position',[.025 .43  .32  .03], ...
          'ToolTipString',...
          ['Any two points whose distance is less than this value are', ...
           ' assumed identical']);

uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{1}, ...
          'Position',[.025 .39  .32  .03], ...
          'ToolTipString',...
          'Minimum constraint violation function value of an infeasible point');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{2}, ...
          'Position',[.025 .36  .32  .03], ...
          'ToolTipString',...
          'Maximum constraint violation function value of any filter point');
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{3}, ...
          'Position',[.025 .33  .32  .03], ...
          'ToolTipString',...
         ['If a discrete neighbor has objective function value within this',...
          ' value of that of the incumbent, then extended polling is', ...
          ' extended polling is performed around this neighbor point']);
uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.other{4}, ...
          'Position',[.025 .30  .32  .03], ...
          'ToolTipString',...
         ['If a discrete neighbor has constraint violation function value', ...
          ' within this value of that of the incumbent, then extended' ...
          ' polling is performed around this neighbor point']);
gui.RSLabel(1) = uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.RS{1}, ...
          'Position',[.025 .26  .32  .03], ...
          'ToolTipString',...
          'Ranking and Selection Parameters');
gui.RSLabel(2) = uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.RS{2}, ...
          'Position',[.025 .23  .32  .03], ...
          'ToolTipString',...
          'Ranking and Selection Parameters');
gui.RSLabel(3) = uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.RS{3}, ...
          'Position',[.025 .20  .32  .03], ...
          'ToolTipString',...
         'Ranking and Selection Parameters');
gui.RSLabel(4) = uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.RS{4}, ...
          'Position',[.025 .17  .32  .03], ...
          'ToolTipString',...
         'Ranking and Selection Parameters');
gui.RSLabel(5) = uicontrol(gui.fig, 'String',gui_var.Labels.Parameters.RS{5}, ...
          'Position',[.025 .14  .32  .03], ...
          'ToolTipString',...
         'Ranking and Selection Parameters');

% Display of NOMADm object values (Problem Name, Parameter Settings, etc.)
gui.problem        = uicontrol(gui.fig, 'String',          '', ...
                                        'FontWeight',      'bold', ...
                                        'ForegroundColor', 'red', ...
                                        'BackgroundColor', [.8 .8 .8], ...
                                        'Position', [.225 .955 .19 .035]);
gui.searchLabel(1) = uicontrol(gui.fig, 'Position', [.14  .87 .27 .03], ...
                                        'ForegroundColor','blue');
gui.searchLabel(2) = uicontrol(gui.fig, 'Position', [.14  .84 .27 .03], ...
                                        'ForegroundColor','blue');
gui.pollStrategy   = uicontrol(gui.fig, 'Position', [.14  .81 .27 .03], ...
                                        'ForegroundColor','blue');
gui.pollOrder      = uicontrol(gui.fig, 'Position', [.14  .78 .27 .03], ...
                                        'ForegroundColor','blue');
gui.pollCenter     = uicontrol(gui.fig, 'Position', [.14  .75 .27 .03], ...
                                        'ForegroundColor','blue');
gui.Term.delta     = uicontrol(gui.fig, 'Position', [.35  .71 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.Term.nIter     = uicontrol(gui.fig, 'Position', [.35  .68 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.Term.nFunc     = uicontrol(gui.fig, 'Position', [.35  .65 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.Term.time      = uicontrol(gui.fig, 'Position', [.35  .62 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.Term.nFails    = uicontrol(gui.fig, 'Position', [.35  .59 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.delta0         = uicontrol(gui.fig, 'Position', [.35  .55 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.deltaMax       = uicontrol(gui.fig, 'Position', [.35  .52 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.meshRefine     = uicontrol(gui.fig, 'Position', [.35  .49 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.meshCoarsen    = uicontrol(gui.fig, 'Position', [.35  .46 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.tolCache       = uicontrol(gui.fig, 'Position', [.35  .43 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.hmin           = uicontrol(gui.fig, 'Position', [.35  .39 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.hmax           = uicontrol(gui.fig, 'Position', [.35  .36 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.ePollXiF       = uicontrol(gui.fig, 'Position', [.35  .33 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.ePollXiH       = uicontrol(gui.fig, 'Position', [.35  .30 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.RS.s0          = uicontrol(gui.fig, 'Position', [.35  .26 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.RS.alpha_const = uicontrol(gui.fig, 'Position', [.35  .23 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.RS.iz_const    = uicontrol(gui.fig, 'Position', [.35  .20 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.RS.alpha_rho   = uicontrol(gui.fig, 'Position', [.35  .17 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');
gui.RS.iz_rho      = uicontrol(gui.fig, 'Position', [.35  .14 .06 .03], ...
                                        'HorizontalAlignment','center', ...
                                        'ForegroundColor','blue');

% Additional ToolTips for NOMADm objects
set([gui.searchLabel(1);gui.searchLabel(2)], ...
    'ToolTipString','The Search screen is accessed in the MADS Menu');
set([gui.pollStrategy;gui.pollOrder;gui.pollCenter], ...
    'ToolTipString','Changing these values is done in the MADS Menu');
set(gui.hmin,'ToolTipString','This value should be kept very small');
set(gui.hmax,'ToolTipString', ...
    'This value must be greater than Min Filter Constraint Violation');
set(gui.TermFlag.delta,'ToolTipString',...
    ['MADS terminates as soon as the mesh size parameter value falls', ...
     ' below this value']);
set([gui.TermFlag.nIter;gui.TermFlag.nFunc; ...
     gui.TermFlag.time; gui.TermFlag.nFails], ...
    'ToolTipString', ...
    'Check the box to activate the corresponding termination criterion');
set(gui.Term.delta,'ToolTipString', ...
    ['This value must be positive and should be kept small as it measures,', ...
     ' the accuracy of the solution']);
set([gui.Term.nIter;gui.Term.nFunc;gui.Term.time;gui.Term.nFails], ...
    'ToolTipString', ...
    'Setting this value to ''Inf'' is the same as unchecking the box');

% Display of 2 plot axes and Stop/Resume Run pushbuttons
gui.axesHistory = axes('Parent',gui.fig, ...
   'Position',            [.52 .605 .43 .33], ...
   'FontUnits',           'normalized', ...
   'FontSize',            .07, ...
   'Box',                 'on', ...
   'Visible',             'off', ...
   'ButtonDownFcn',       gui_func.copyPlot);
gui.axesFilter = axes('Parent',gui.fig, ...
   'Position',            [.62 .13 .33 .33], ...
   'FontUnits',           'normalized', ...
   'FontSize',            .07, ...
   'Box',                 'on', ...
   'Visible',             'off', ...
   'ButtonDownFcn',       gui_func.copyPlot);
gui.stopRun = uicontrol(gui.fig, ...
   'Style',               'pushbutton', ...
   'String',              'Stop Run', ...
   'Tag',                 'StopRun', ...
   'FontWeight',          'bold', ...
   'HorizontalAlignment', 'center', ...
   'Position',            [.42 .50 .11 .045], ...
   'FontSize',            .5, ...
   'BackgroundColor',     [.8 .8 .8], ...
   'Visible',             'off', ...
   'Interruptible',       'off', ...
   'UserData',            0, ...
   'Callback',            'set(gcbo,''UserData'',1);');
gui.resumeRun = uicontrol(gui.fig, ...
   'Style',               'pushbutton', ...
   'String',              'Resume Run', ...
   'FontWeight',          'bold', ...
   'HorizontalAlignment', 'center', ...
   'Position',            [.42 .455 .11 .045], ...
   'FontSize',            .5, ...
   'BackgroundColor',     [.8 .8 .8], ...
   'Visible',             'off', ...
   'Interruptible',       'off', ...
   'Callback',            {gui_func.runMADS, 2});

% Set up Main Menu
gui.Menu.Problem = uimenu(gui.fig,'Position',1,'Label','Problem');
gui.Menu.MADS    = uimenu(gui.fig,'Position',2,'Label','MADS');
gui.Menu.Options = uimenu(gui.fig,'Position',3,'Label','Options');
gui.Menu.Session = uimenu(gui.fig,'Position',4,'Label','Session');
gui.Menu.Cache   = uimenu(gui.fig,'Position',5,'Label','Cache');
gui.Menu.Run     = uimenu(gui.fig,'Position',6,'Label','Run');
gui.Menu.Results = uimenu(gui.fig,'Position',7,'Label','Results');
gui.Menu.Help    = uimenu(gui.fig,'Position',8,'Label','Help');

% Define Problem Menu
gui.ProblemMenu.new     = uimenu(gui.Menu.Problem, 'Position',  1, ...
   'Label',       'New Optimization Problem', ...
   'Accelerator', 'N', ...
   'Enable',      'on', ...
   'Callback',    gui_func.selectProblem);
for k = 1:length(gui_var.Types.file)
    gui.ProblemMenu.edit(k) = uimenu(gui.Menu.Problem, 'Position', k+1, ...
      'Label',    ['Edit ', gui_var.Labels.file{k}, ' File'], ...
      'Visible',  'off', ...
      'CallBack',  gui_func.editFile);
end
gui.ProblemMenu.quit = uimenu(gui.Menu.Problem, ...
   'Position',    length(gui_var.Types.file)+2, ...
   'Label',       'Quit NOMADm', ...
   'Accelerator', 'Q', ...
   'Callback',    {gui_func.quitGUI, gui_var});
set([gui.ProblemMenu.edit(1);gui.ProblemMenu.quit],'Separator','on');

% Define MADS Menu
gui.MADSMenu.editTermParam  = uimenu(gui.Menu.MADS, 'Position', 1, ...
   'Label',       'Edit Termination Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.term, ...
                  [gui.Term.delta;gui.Term.nIter;gui.Term.nFunc; ...
                  gui.Term.time;gui.Term.nFails]}, ...
  'Callback',     gui_func.editParameters);
gui.MADSMenu.editMeshParam  = uimenu(gui.Menu.MADS, 'Position', 2, ...
   'Label',       'Edit Mesh Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.mesh, ...
                  [gui.delta0;gui.deltaMax;gui.meshRefine; ...
                  gui.meshCoarsen;gui.tolCache]}, ...
   'Callback',    gui_func.editParameters);
gui.MADSMenu.editOtherParam = uimenu(gui.Menu.MADS, 'Position', 3, ...
   'Label',       'Edit Filter/MVP Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.other, ...
                  [gui.hmin; gui.hmax; gui.ePollXiF; gui.ePollXiH]}, ...
   'Callback',    gui_func.editParameters);
gui.MADSMenu.editRSParam    = uimenu(gui.Menu.MADS, 'Position', 4, ...
   'Label',       'Edit Ranking & Selection Parameters', ...
   'UserData',    {gui_var.Labels.Parameters.RS, ...
                  [gui.RS.s0; gui.RS.alpha_const; gui.RS.iz_const; ...
                   gui.RS.alpha_rho; gui.RS.iz_rho]}, ...
   'Callback',    gui_func.editParameters);
   
gui.MADSMenu.editSearch     = uimenu(gui.Menu.MADS, 'Position', 5, ...
   'Separator',   'on', ...
   'Label',       'Select Search Strategies', ...
   'Callback',    gui_func.search_gui);
gui.MADSMenu.pollStrategies = uimenu(gui.Menu.MADS, 'Position', 6, ...
   'Separator',   'on', ...
   'Label',       'Select Poll Directions');
gui.MADSMenu.pollOrders     = uimenu(gui.Menu.MADS, 'Position', 7, ...
   'Label',       'Select Poll Order');
gui.MADSMenu.pollCenters    = uimenu(gui.Menu.MADS, 'Position', 8, ...
   'Label',       'Select Poll Center');
gui.MADSMenu.pollComplete   = uimenu(gui.Menu.MADS, 'Position', 9, ...
   'Label',       'Complete Polling', ...
   'Callback',    'umtoggle(gcbo);');

% Define MADS Submenus
for k = 1:length(gui_var.Labels.pollStrategy)
   gui.MADSMenu.pollStrategy(k) = uimenu(gui.MADSMenu.pollStrategies, ...
      'Position', k, ...
      'Label',    gui_var.Labels.pollStrategy{k}, ...
      'UserData', gui.pollStrategy, ...
      'Callback', {gui_func.selectPoll,'Strategy'});
   if strncmp(gui_var.Types.poll{k}, 'Gradient', 8)
      set(gui.MADSMenu.pollStrategy(k),'Enable','off');
   end
end
for k = 1:length(gui_var.Labels.pollCenter)
   gui.MADSMenu.pollCenter(k) = uimenu(gui.MADSMenu.pollCenters, ...
      'Position', k, ...
      'Label',    gui_var.Labels.pollCenter{k}, ...
      'UserData', gui.pollCenter, ...
      'Callback', {gui_func.selectPoll,'Center'});
end
for k = 1:length(gui_var.Labels.pollOrder)
   gui.MADSMenu.pollOrder(k) = uimenu(gui.MADSMenu.pollOrders, ...
      'Position', k, ...
      'Label',    gui_var.Labels.pollOrder{k}, ...
      'UserData', gui.pollOrder, ...
      'Callback', {gui_func.selectPoll,'Order'});
end

% Define Options Menu
gui.OptionsMenu.useFilter         = uimenu(gui.Menu.Options, 'Position',1, ...
   'Label',       'Use Filter for Nonlinear Constraints', ...
   'Callback',    'umtoggle(gcbo);');
gui.OptionsMenu.removeRedundancy  = uimenu(gui.Menu.Options, 'Position',2, ...
   'Label',       'Discard Redundant Linear Constraints', ...
   'Callback',    'umtoggle(gcbo);');
gui.OptionsMenu.runStochastic     = uimenu(gui.Menu.Options, 'Position',3, ...
   'Label',       'Run as Stochastic Optimization Problem', ...
   'Callback',    gui_func.toggleStochastic);
gui.OptionsMenu.scaleMenu         = uimenu(gui.Menu.Options, 'Position',4, ...
   'Label',       'Scaling of Mesh Directions');
gui.OptionsMenu.accelerate        = uimenu(gui.Menu.Options, 'Position',5, ...
   'Label',       'Accelerate Convergence', ...
   'Callback',    'umtoggle(gcbo);');
gui.OptionsMenu.TermFlag.relative = uimenu(gui.Menu.Options, 'Position',6, ...
   'Label',       'Use Relative Termination Tolerance', ...
   'Callback',    'umtoggle(gcbo);');
gui.OptionsMenu.plotHistory1      = uimenu(gui.Menu.Options, 'Position',7, ...
   'Label',       'Plot History', 'Separator',   'on', ...
   'Callback',    'umtoggle(gcbo);');
gui.OptionsMenu.plotHistory2      = uimenu(gui.Menu.Options, 'Position',8, ...
   'Label',       'Plot History (Real-time)', ...
   'Callback',    'umtoggle(gcbo);');
gui.OptionsMenu.plotFilter        = uimenu(gui.Menu.Options, 'Position',9, ...
   'Label',       'Plot Filter (Real-time)', ...
   'Callback',    'umtoggle(gcbo);');

% Define Scale SubMenu of the Options Menu choice Scale
for k = 1:length(gui_var.Labels.scale)
   gui.scaleMenu(k) = uimenu(gui.OptionsMenu.scaleMenu, 'Position', k, ...
      'Label',    gui_var.Labels.scale{k}, ...
      'Callback', ['set(get(get(gcbo,''Parent''),''Children''),', ...
                   ' ''Checked'',''off''); umtoggle(gcbo);']);
end

% Define Session Menu
gui.SessionMenu.load   = uimenu(gui.Menu.Session, 'Position', 1, ...
   'Label',       'Load Options from Session File', ...
   'Enable',      'off', ...
   'Callback',    gui_func.loadSession);
gui.SessionMenu.save   = uimenu(gui.Menu.Session, 'Position', 2, ...
   'Label',       'Save Options to Session File', ...
   'Enable',      'off', ...
   'Callback',    gui_func.saveSession);
gui.SessionMenu.reset  = uimenu(gui.Menu.Session, 'Position', 3, ...
   'Label',       'Reset Options to Defaults', ...
   'Callback',    gui_func.resetSession);
gui.SessionMenu.delete = uimenu(gui.Menu.Session, 'Position', 4, ...
   'Label',       'Delete Session File', ...
   'Enable',      'off', ...
   'UserData',   {'Session File',gui_var.FileExt.S,gui.SessionMenu.load},...
   'Callback',    gui_func.deleteFile);

% Define Cache Menu
gui.CacheMenu.load   = uimenu(gui.Menu.Cache, 'Position', 1, ...
   'Label',       'Use Pre-Existing Cache File', ...
   'Enable',      'off', ...
   'Callback',    'umtoggle(gcbo);');
gui.CacheMenu.count  = uimenu(gui.Menu.Cache, 'Position', 2, ...
   'Label',       'Count Cache Points as Function Calls', ...
   'Enable',      'off', ...
   'Callback',    'umtoggle(gcbo);');
gui.CacheMenu.save   = uimenu(gui.Menu.Cache, 'Position', 3, ...
   'Label',       'Save Run Results to Cache File', ...
   'Accelerator', 'S', ...
   'Enable',      'off', ...
   'Callback',    gui_func.saveCache);
gui.CacheMenu.delete = uimenu(gui.Menu.Cache, 'Position', 4, ...
   'Label',       'Delete Pre-Existing Cache File', ...
   'Enable',      'off', ...
   'UserData',    {'Cache File', gui_var.FileExt.C, ...
                  [gui.CacheMenu.load; gui.CacheMenu.count]}, ...
   'Callback',    gui_func.deleteFile);

% Define Run Menu
gui.RunMenu.exec         = uimenu(gui.Menu.Run,'Position',1, ...
   'Label',       'Execute Next Run', ...
   'Accelerator', 'R', ...
   'Enable',      'off', ...
   'Callback',    {gui_func.runMADS, 0});
gui.RunMenu.resume       = uimenu(gui.Menu.Run,'Position',2, ...
   'Label',       'Resume a Stopped Run', ...
   'Enable',      'off', ...
   'Callback',    {gui_func.runMADS, 2});
gui.RunMenu.restart  = uimenu(gui.Menu.Run,'Position',3, ...
   'Label',       'Restart Run from Current Point', ...
   'Enable',      'off', ...
   'Callback',    {gui_func.runMADS, 1});
gui.RunMenu.oneIteration = uimenu(gui.Menu.Run,'Position',4, ...
   'Separator',   'on', ...
   'Label',       'Run One Iteration at a Time', ...
   'Checked',     'off', ...
   'Enable',      'off', ...
   'Callback',    'umtoggle(gcbo);');
gui.RunMenu.execFeasible = uimenu(gui.Menu.Run,'Position',5, ...
   'Label',       'Run Only Until Feasible', ...
   'Checked',     'off', ...
   'Enable',      'off', ...
   'Callback',    'umtoggle(gcbo);');
gui.RunMenu.clear        = uimenu(gui.Menu.Run,'Position',6, ...
   'Separator',   'on', ...
   'Label',       'Clear Previous Runs', ...
   'Enable',      'off', ...
   'Callback',    gui_func.clearRuns);

% Define Results Menu
set(gui.Menu.Results,'Visible','off');
for k = 1:gui_var.runMax
   gui.ResultsMenu(k) = uimenu(gui.Menu.Results,'Position',k, ...
      'Label',          ['View Run # ' int2str(k)], ...
      'Visible',        'off', ...
      'ForegroundColor',gui_var.Types.plotColors(k), ...
      'Accelerator',    int2str(mod(k,gui_var.runMax)), ...
      'Callback',       gui_func.results_gui);
end

% Define Help Menu
uimenu(gui.Menu.Help,'Position',1, ...
   'Label',      'NOMADm Help', ...
   'Accelerator','H', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.nomadm,'file')}, ...
   'UserData',   gui_var.HelpDoc.nomadm,...
   'Callback',   gui_func.getHelp);
uimenu(gui.Menu.Help,'Position',2, ...
   'Label',      'DACE Toolbox Help', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.dace,'file')}, ...
   'UserData',   gui_var.HelpDoc.dace,...
   'Callback',   gui_func.getHelp);
uimenu(gui.Menu.Help,'Position',3, ...
   'Label',      'N-W Toolbox Help', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.nw,'file')}, ...
   'UserData',   gui_var.HelpDoc.nw,...
   'Callback',   gui_func.getHelp);
uimenu(gui.Menu.Help,'Position',4, ...
   'Label',      'CMA-ES Toolbox Help', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.cmaes,'file')}, ...
   'UserData',   gui_var.HelpDoc.cmaes,...
   'Callback',   gui_func.getHelp);
uimenu(gui.Menu.Help,'Position',5, ...
   'Separator',  'on', ...
   'Label',      'View List of Version Changes', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.changes,'file')}, ...
   'UserData',   gui_var.HelpDoc.changes,...
   'Callback',   gui_func.getHelp);
uimenu(gui.Menu.Help,'Position',6, ...
   'Label',      'View GNU Public License', ...
   'Visible',    onoff{2-~~exist(gui_var.HelpDoc.license,'file')}, ...
   'UserData',   gui_var.HelpDoc.license,...
   'Callback',   gui_func.getHelp);
uimenu(gui.Menu.Help,'Position',7, ...
   'Separator',  'on', ...
   'Label',      'About NOMADm', ...
   'Callback',   {gui_func.getAbout,gui_var.version,gui_var.Labels.coolCats});
return

%*******************************************************************************
% BEGINNING OF SIMPLE NOMADm GUI CALLBACK FUNCTIONS
%*******************************************************************************

%*******************************************************************************
% selectProblem: Choose an optimization problem.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Problem-->New Optimization Problem)
% Calls:     clearRuns, resetSession
% VARIABLES:
%  pFile             = filename of current optimization problem
%  pPath             = pathname of current optimization problem
%  gui_var           = structure of all GUI variables
%    .path           =   path of current optimization problem
%    .problemExt     =   filename extension of optimization problem 
%    .problem        =   name of current optimization problem
%    .FileExt        =   filename suffixes
%      .O            =     Omega file suffix
%      .I            =     initial points file suffix
%      .N            =     discrete neighbors file suffix
%      .P            =     user parameter file suffix
%      .C            =     Cache file suffix
%      .S            =     Session file suffix
%    .noGrad         =   flag indicating initial availability of gradients
%  gui               = structure of all GUI object handles
%    .runStatus      =   handle for text bar at bottom of figure window
%    .RunMenu        =   handles for Run menu items
%    .ProblemMenu    =   handles for Problem menu items
%    .MADSMenu       =   handles for MADS menu items
%    .SessionMenu    =   handles for Session menu items
%    .CacheMenu      =   handles for CacheMenu items
%  nameProblem       = name of current optimization problem
%  isConstrained     = flags problem as having linear constraints
%  hasInitGuess      = flags problem as having an initial guess
%  hasUserParam      = flags problem as having user parameters
%  isMVP             = flags problem as being an MVP
%  existCache        = flags problem as having an existing Cache file
%  existSession      = flags problem as having an existing Session file
%  pLabels           = strings of Poll strategy labels
%*******************************************************************************
function selectProblem(h,event)

[pFile,pPath] = uigetfile({'*.m',        'Matlab M-files (*.m)';      ...
                           '*.f; *.for', 'Fortran files (*.f,*.for)'; ...
                           '*.c; *.C',   'C/C++ files (*.c,*.C)'},    ...
                           'Select an Optimization Problem file');
if (pFile)

   % Reset all variables and retrieve GUI application data
   clearRuns(1,0);
   resetSession(1,0);
   gui_var = getappdata(gcbf,'gui_var');
   gui     = getappdata(gcbf,'gui');

   % Set the path, removing only directories not there at startup
   if gui_var.newPath
      rmpath(gui_var.path);
   end
   gui_var.newPath = isempty(strfind(upper(path),upper(pPath)));
   gui_var.path    = pPath;
   if gui_var.newPath
      addpath(gui_var.path);
   end

   % Load New Problem
   [nameProblem, gui_var.problemExt] = strtok(pFile,'.');
   set(gui.problem,             'String',nameProblem);
   set(gui.runStatus,           'String','No runs performed');
   set(gui.RunMenu.exec,        'Enable','on');
   set(gui.RunMenu.oneIteration,'Enable','on');
   set(gui.RunMenu.execFeasible,'Enable','on');

   % Allow editing only of files that exist
   isConstrained = exist([nameProblem,gui_var.FileExt.O,'.m'],'file');
   isMVP         = exist([nameProblem,gui_var.FileExt.N,'.m'],'file');
   hasInitGuess  = exist([nameProblem,gui_var.FileExt.I,'.m'],'file');
   hasUserParam  = exist([nameProblem,gui_var.FileExt.P,'.m'],'file');
   existCache    = exist([nameProblem,gui_var.FileExt.C],'file');
   existSession  = exist([nameProblem,gui_var.FileExt.S],'file');
   onoff = {'on','off'};
   set(gui.ProblemMenu.edit(1),    'Visible', 'on');
   set(gui.ProblemMenu.edit(2),    'Visible', onoff{1+~hasInitGuess});
   set(gui.ProblemMenu.edit(3),    'Visible', onoff{1+~isConstrained});
   set(gui.ProblemMenu.edit(4),    'Visible', onoff{1+~isMVP});
   set(gui.ProblemMenu.edit(5),    'Visible', onoff{1+~hasUserParam});
   set(gui.MADSMenu.pollOrder(end),'Enable',  onoff{1+~hasUserParam});
   set(gui.SessionMenu.load,       'Enable',  onoff{1+~existSession});
   set(gui.SessionMenu.save,       'Enable',  'on');
   set(gui.SessionMenu.delete,     'Enable',  onoff{1+~existSession});
   set(gui.CacheMenu.load,         'Enable',  onoff{1+~existCache});
   set(gui.CacheMenu.count,        'Enable',  onoff{1+~existCache});
   set(gui.CacheMenu.delete,       'Enable',  onoff{1+~existCache});      
      
   % Allow gradient options if functions file returns enough arguments
   try
      gui_var.noGrad = (abs(nargout(nameProblem)) <= 2);
   catch
      gui_var.noGrad = 0;
   end
   pLabels = char(gui_var.Labels.pollStrategy);
   set(gui.MADSMenu.pollStrategy(find(pLabels(:,1) == 'G')), ...
                               'Enable',onoff{1+gui_var.noGrad});
   setappdata(gcbf,'gui_var',gui_var);
end
return

%*******************************************************************************
% editFile: Edit an optimization problem file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Problem-->Edit XXXXX File)
% VARIABLES:
%  k             = Problem Menu position number of selected file
%  ext           = filename extension of file to be edited
%  filename      = full file name of file to be edited
%  gui_var       = structure of all GUI variables
%    .path       =   path of current optimization problem
%    .problemExt =   filename extension of optimization problem 
%    .Types.file =   strings of file suffixes
%  gui.problem   = name of current optimization problem
%*******************************************************************************
function editFile(h,event)

[h,fig] = gcbo;
gui_var = getappdata(fig,'gui_var');
gui     = getappdata(fig,'gui');

k = get(h,'Position') - 1;
if (k == 1)
   ext = gui_var.problemExt;
else
   ext = '.m';
end
filename = fullfile(gui_var.path, [get(gui.problem,'String'),...
                    gui_var.Types.file{k},ext]);
edit(filename);
clear(filename);
return

%*******************************************************************************
% quitGUI: End program session.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Problem-->Quit NOMADm)
%            nomadm_gui (Close Request Function for figure window)
% VARIABLES:
%  gui_var    = structure of all NOMADm variables
%    .newPath = logical showing if the problem path is on the Matlab path
%    .path    = path of current optimization problem
%*******************************************************************************
function quitGUI(h, event, gui_var)

if gui_var.newPath
   rmpath(gui_var.path);
end
if ishandle(gcbf)
   delete(gcbf);
end
clear variables;
clear functions;
return

%*******************************************************************************
% editParameters: Edit parameters.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Edit Mesh Parameters)
%            nomadm_gui (callback: MADS-->Edit Filter/MVP Parameters)
%            nomadm_gui (callback: MADS-->Edit Termination Parameters)
% VARIABLES:
%  param   = User data storage
%  header  = title of input dialog box
%  labels  = text labels in input dialog box
%  handles = handles to appropriate GUI figure objects
%  defAns  = default answers that appear in input dialog box
%  newAns  = new answers that appear in input dialog box
%*******************************************************************************
function editParameters(h,event)

header  = get(gcbo,'Label');
param   = get(gcbo,'UserData');
[labels, handles] = deal(param{:});
defAns  = get(handles, 'String');
newAns  = inputdlg(labels,header,ones(length(defAns),1)*[1,45],defAns);
if (~isempty(newAns))
   set(handles, {'String'}, newAns);
end
return

%*******************************************************************************
% selectPoll: Select a Poll Strategy, Poll Center, or Poll Order.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Select Poll Directions), 
%            nomadm_gui (callback: MADS-->Select Poll Center),
%            nomadm_gui (callback: MADS-->Select Poll Order), 
% VARIABLES:
%  type            = string containing "Strategy", "Center", or "Order"
%  gui_var.Choice  = user choices
%  gui             = structure of NOMADm object handles
%*******************************************************************************
function selectPoll(h,event,type)

% Update Poll option and screen
gui_var = getappdata(gcbf,'gui_var');
gui_var.Choice.(['poll',type]) = get(gcbo,'Position');
set(get(gcbo,'UserData'),'String',get(gcbo,'Label'));
setappdata(gcbf,'gui_var',gui_var);
return

%*******************************************************************************
% toggleStochastic: Turn on and off flag for Running as a stochastic problem.
% ------------------------------------------------------------------------------
% Called by: updateScreen, nomadm_gui (callback: MADS-->
%                                      Run as Stochastic Optimization Problem)
% VARIABLES:
%  gui                          = structure of NOMADm object handles
%    .OptionsMenu.runStochastic = menu option for running stochastic problem
%    .MADSMenu.editRSParam      = menu option for editing R&S parameters
%    .RSPanel                   = handle of the R&S parameter display panel
%    .RS                        = structure of handles for R&S parameter values
%    .RSLabel                   = handles of R&S parameter display labels
%  z                            = value of the menu item toggle switch
%  flag                         = conversion of z to "on" or "off"
%*******************************************************************************
function toggleStochastic(h,event)

% Determine if stochastic option is turned on or off
gui     = getappdata(gcf,'gui');
gui_var = getappdata(gcf,'gui_var');
z = umtoggle(gui.OptionsMenu.runStochastic);
if z == 1
   flag = 'on';
else
   flag = 'off';
end;

% Turn on/off R&S displays
set(gui.MADSMenu.editRSParam,   'Enable', flag);
set([gui.RSPanel;gui.RSLabel'], 'Visible', flag);
if ~gui_var.newMatlab
   set([gui.RSBorder], 'Visible', flag);
end
rsfield = fieldnames(gui.RS);
for k = 1:length(rsfield)
   set(gui.RS.(rsfield{k}), 'Visible', flag);
end
return

%*******************************************************************************
% loadSession: Load session options from previously saved file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Load Options from Session File)
% Calls:     updateScreen
% VARIABLES:
%  sessionFile     = full path and name of session file
%  gui_var         = structure of NOMADm variables
%    .path         = path of current optimization problem
%    .FileExt.S    = filename extension of session file
%    .Choice       = current user choices
%    .Options      = current option settings
%  gui.problem     = object handle for Optimization Problem Name
%  Session.Choice  = previously saved user choices 
%  Session.Options = previously saved option settings
%*******************************************************************************
function loadSession(h,event)

gui_var = getappdata(gcf,'gui_var');
gui     = getappdata(gcf,'gui');

sessionFile = fullfile(gui_var.path, ...
              [get(gui.problem,'String'), gui_var.FileExt.S]);
if exist(sessionFile,'file')
   load(sessionFile);
   gui_var.Choice  = Session.Choice;
   gui_var.Options = Session.Options;
   updateScreen(gui_var.Choice,gui_var.Options);
   setappdata(gcf,'gui_var',gui_var);
end
return

%*******************************************************************************
% resetSession: Reset MADS Parameters to Default values.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Reset Options to Defaults),
%            nomadm
% Calls:     updateScreen
% VARIABLES:
%  gui_var.Choice           = current user choices
%  gui_var.Options          = current option settings
%  gui_var.Defaults.Choice  = default user choice settings
%  gui_var.Defaults.Options = default options settings
%*******************************************************************************
function resetSession(h,event)

gui_var = getappdata(gcf,'gui_var');
gui_var.Choice  = gui_var.Defaults.Choice;
gui_var.Options = gui_var.Defaults.Options;
updateScreen(gui_var.Choice,gui_var.Options);
setappdata(gcf,'gui_var',gui_var);
return

%*******************************************************************************
% saveSession: Save session options to file for future retrieval.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Save Options to Session File)
% Calls:     loadMADS
% VARIABLES:
%  sessionFile            = full path and name of session file
%  gui_var                = structure of all GUI variables
%    .path                =   path of current optimization problem
%    .FileExt.S           =   filename extension of session file
%    .Choice              =   current user choices
%    .Options             =   current option settings
%  gui                    = structure of GUI object handles
%    .problem             =   name of current optimization problem
%    .SessionMenu         =   handles for Session menu items
%    .OptionsMenu         =   handles for Options menu items
%    .runStatus           =   handle for figure window status bar
%  Session                = previously saved options and parameters
%    .Problem             =   optimization problem data
%    .Choice              =   user choices 
%    .Options             =   option settings
%      .Term              =     termination criteria
%        .iter            =       number of iterations
%        .func            =       number of function evaluations
%        .time            =       CPU time
%        .fails           =       number of consecutive Poll failures
%      .TermFlag.relative =     termination relative to initial mesh size
%*******************************************************************************
function saveSession(h,event)

gui     = getappdata(gcbf,'gui');
gui_var = getappdata(gcbf,'gui_var');

sessionFile = [get(gui.problem,'String'), gui_var.FileExt.S];
[Session.Problem,Session.Options] = loadMADS(gui,gui_var);
Session.Options.Term.nIter  = get(gui.Term.nIter,  'String');
Session.Options.Term.nFunc  = get(gui.Term.nFunc,  'String');
Session.Options.Term.time   = get(gui.Term.time,   'String');
Session.Options.Term.nFails = get(gui.Term.nFails, 'String');
umtoggle(gui.OptionsMenu.TermFlag.relative);
Session.Options.TermFlag.relative = umtoggle(gui.OptionsMenu.TermFlag.relative);
Session.Choice = gui_var.Choice;
save([gui_var.path, sessionFile],'Session');
set(gui.runStatus,'String', ...
   ['Session Options saved to file, ',sessionFile]);
set([gui.SessionMenu.load; gui.SessionMenu.delete],'Enable','on');
return

%*******************************************************************************
% saveCache: Save MADS run results in a .mat file for later retrieval.
% ------------------------------------------------------------------------------
% -
% Called by: nomadm_gui (callback: Cache-->Save Run Results to Cache File),
%            loadMADS
% VARIABLES:
%  gui             = structure of all GUI object handles
%    .fig          =   handle for the GUI figure window
%    .problem      =   name of current optimization problem
%    .runStatus    =   handle for the GUI figure window status bar
%    .CacheMenu    =   handles for Cache menu items
%  gui_var         = structure of all GUI variables
%    .runCount     =   MADS run number
%    .path         =   path of current optimization problem
%    .FileExt.C    =   default Cache filename suffix 
%  CName           = name of Cache file
%  RunSet(1).Cache = structure of MADS run data
%  Cache           = storage of Cache for use by MADS
%*******************************************************************************
function saveCache(h,event)

if isappdata(gcbf,'RunSet')
   gui       = getappdata(gcbf,'gui');
   gui_var   = getappdata(gcbf,'gui_var');
   RunSet    = getappdata(gcbf,'RunSet');
   setptr(gcf,'watch');
   set(gui.runStatus, 'String', ...
       ['Run # ',int2str(gui_var.runCount),' Cache being saved']);
   CName = [get(gui.problem, 'String'), gui_var.FileExt.C];
   Cache = RunSet(1).Cache;
   save([gui_var.path, CName],'Cache');
   set([gui.CacheMenu.load;gui.CacheMenu.count;gui.CacheMenu.delete], ...
       'Enable','on');
   set(gui.runStatus, 'String', ...
       ['Run # ',int2str(gui_var.runCount),' Cache saved in ',CName]);
   setptr(gcbf,'arrow');
else
   error('Cannot save.  Cache not found');
end
return

%*******************************************************************************
% deleteFile: Delete a Session or Cache file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Session-->Delete Session File), 
%            nomadm_gui (callback: Cache-->Delete Pre-Existing Cache File)
% VARIABLES:
%  param          = user data that identifies file to be deleted
%  fileID         = string identifing file to be deleted
%  fileExt        = filename suffix of file to be deleted
%  handles        = handles of menu choices to be disabled
%  file           = full path and name of file to be deleted
%  deleteFile     = logical for deleting file
%  gui_var.path   = path of current optimization problem
%  gui            = structure of all GUI object handles
%    .runStatus   =   object handle for figure window status bar
%    .problem     =   name of current optimization problem
%*******************************************************************************
function deleteFile(h,event)

gui_var = getappdata(gcbf,'gui_var');
gui     = getappdata(gcbf,'gui');
param   = get(gcbo,'UserData');
[fileID,fileExt,handles] = deal(param{:});
file  = fullfile(gui_var.path,[get(gui.problem,'String'),fileExt]);
deleteFile = questdlg(['Are you sure you want to delete ',file,'?'], ...
                      ['Delete ',fileID,'?'],'No');
if strcmp(deleteFile, 'Yes')
   delete(file);
   set(gui.runStatus,'String', [fileID,', ',file,', has been deleted']);
   set([handles; gcbo],'Enable','off');
end
return

%*******************************************************************************
% clearRuns: Clear figure window and reset all the appropriate variables.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Run-->Clear Previous Runs),
%            selectProblem
% VARIABLES:
%  gui               = structure of all GUI object handles
%    .fig            =   handle for GUI figure window
%    .stopRun        =   handle for Stop Run pushbutton
%    .resumeRun      =   handle for Resume Run pushbutton
%    .axesHistory    =   handle for History plot
%    .axesFilter     =   handle for Filter plot
%    .ResultsMenu    =   handles for Results menu items
%    .RunMenu        =   handles for Run menu items
%    .CacheMenu      =   handles for Cache menu items
%    .runStatus      =   handle for figure window status bar
%  gui_var           = structure of all GUI variables
%    .runMax         =   maximum nmber of MADS runs
%    .runCount       =   MADS run counter
%  RunSet            = structure of MADS run data
%*******************************************************************************
function clearRuns(h,event)

% Delete application data
if isappdata(gcbf,'RunSet')
   rmappdata(gcbf,'RunSet');
end
if isappdata(0,'CACHE')
   rmappdata(0,'CACHE');
end
if isappdata(0,'sCACHE')
   rmappdata(0,'sCACHE');
end
if isappdata(0,'surrogate')
   rmappdata(0,'surrogate');
end

% Reset properties of gui object handles
gui = getappdata(gcbf,'gui');
set([gui.RunMenu.clear; gui.CacheMenu.save],             'Enable', 'off');
set([gui.RunMenu.resume;gui.RunMenu.restart],            'Enable', 'off');      
set([gui.stopRun;       gui.resumeRun],                  'Visible','off');      
set([gui.Menu.Results;  gui.ResultsMenu'],               'Visible','off');
set([gui.axesFilter;    get(gui.axesFilter, 'Children')],'Visible','off');
set([gui.axesHistory;   get(gui.axesHistory,'Children')],'Visible','off');
set( gui.axesHistory,   'NextPlot','replace');
set( gui.runStatus,     'String',  'All runs cleared');

% Reset the run counter
gui_var = getappdata(gcbf,'gui_var');
gui_var.runCount = 0;
setappdata(gcbf,'gui_var',gui_var);

% Reset everything else
cla;
setptr(gcbf,'arrow');
clear functions;
return

%*******************************************************************************
% getHelp: View a Help file.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Help-->NOMADm Help),
%            nomadm_gui (callback: Help-->DACE Toolbox Help),
%            nomadm_gui (callback: Help-->NW Toolbox Help),
%            nomadm_gui (callback: Help-->CMA-ES Toolbox Help),
%            nomadm_gui (callback: Help-->View List of Version CHanges),
%            nomadm_gui (callback: Help-->View GNU Public License),
% VARIABLES:
%  helpfile  = name of the help file to be displayed
%  errormsg  = message to display when error is flagged
%  errorflag = error flag for viewing failure
%*******************************************************************************
function getHelp(h,event)

helpfile = which(get(gcbo,'UserData'));
if isempty(helpfile)
   errormsg = ['Error: Help file not found on the Matlab path.  ', ...
               'Use "File-->Set Path" to add location to Matlab path.'];
   errordlg(errormsg);
end
errorflag = web(['file:///', helpfile], '-browser');
if (errorflag), errordlg('Error: Browser or help file not found'); end
return

%*******************************************************************************
% getAbout: View author information.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Help-->About NOMADm)
% VARIABLES:
%    pversion     = current version number of this software
%    contributors = list of people who have contributed to this software
%*******************************************************************************
function getAbout(h,event,pversion,contributors)

msgbox({['NOMADm, version ',pversion], ...
        'Copyright (c) 2001-2005 by Mark A. Abramson','', ...
        'Special Thanks to:','', ...
         contributors{:}},get(gcbo,'Label'));
return

%*******************************************************************************
% copyPlot: Copy history or filter plot to a new screen.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (ButtonDown Function: History and Filter Plots)
% VARIABLES:
%  h1 = handle to new figure window
%  h2 = handle to new axes on figure window
%*******************************************************************************
function copyPlot(h,event)

h1 = figure;
h2 = copyobj(gcbo,h1);
set(h2,'Position','default','ButtonDownFcn','');
return

%*******************************************************************************
% END OF SIMPLE NOMADm GUI CALLBACK FUNCTIONS
%*******************************************************************************

%*******************************************************************************
% runMADS:  Assign GUI input fields to MADS variables and run MADS.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Run-->Execute Next Run)
%            nomadm_gui (callback: Run-->Resume a Stopped Run)
%            nomadm_gui (callback: Run-->Restart Run from Current Point)
% Calls:     loadMADS, saveCache, mads
% VARIABLES:
%  restart           = flag for resuming previous run
%  gui_var           = structure containing all GUI variables
%    .runCount       =   current MADS run number
%    .runMax         =   maximum allowed MADS run number
%    .Options.runUntilFeasible = flag for running until a feasible poi
%  gui_func          = structure of all NOMADm function handles
%  gui               = structure of GUI object handles
%    .fig            =   handle for the main figure window
%    .stopRun        =   handle for Stop Run pushbutton
%    .resumeRun      =   handle for Resume Run pushbutton
%    .axesFilter     =   handle for Filter plot
%    .axesHistory    =   handle for History plot
%    .RunMenu        =   handles for Run menu items
%    .ResultsMenu    =   handles for Results menu items
%    .CacheMenu      =   handles for Cache menu items
%    .runStatus      =   handle for figure window status bar
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
%    .runCount       =   current MADS run number
%  iterate0          = initial iterate
%  BestF             = best feasible iterate found by MADS
%  BestI             = least infeasible iterate found by MADS
%  RunStats          = structure containing MADS run statistics
%*******************************************************************************
function runMADS(h,event,restart)

% Get application data
fig       = findobj(0,'Tag','NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');
gui_func  = getappdata(fig,'gui_func');
gui       = getappdata(fig,'gui');
if isappdata(fig,'RunSet')
   RunSet = getappdata(fig,'RunSet');
end

% Set flag for running only until a feasible point is found
if (restart == -1)
   gui_var.Options.runUntilFeasible = 1;
end

lasterr('NOMADm interrupted by user');
try
catch
end
   % Change figure window, as appropriate
   if (gui_var.runCount >= gui_var.runMax)
      error('Too many MADS runs without clearing.');
   end
   setptr(fig,'watch');
   set(gui.stopRun,   'Visible','on','UserData',0);
   set(gui.resumeRun, 'Visible','off');
   set(gui.RunMenu.clear,'Enable','on');
   set([gui.axesFilter; get(gui.axesFilter,'Children')], 'Visible', 'off');
   set(gui.axesFilter,'ButtonDownFcn','');
   drawnow;

   % Load MADS input data
   [Problem,Options] = loadMADS(gui,gui_var);
   Options.Search(Options.nSearches+1:end) = [];
   
   % Get initial point
   if (restart > 0)
      iterate0 = [RunSet(gui_var.runCount).BestF, ...
                  RunSet(gui_var.runCount).BestI];
   elseif (~restart && exist(Problem.File.I,'file'))
      iterate0 = feval(Problem.File.I);
   else
      error(['Cannot find: ', Problem.File.I, '.m']);
   end

   % Set up and run MADS algorithm and store output
   if (restart == 2)
      saveCache(1,0);
      set(gui.axesHistory, 'NextPlot','replacechildren');
      Options.loadCache  = 1;
      Options.countCache = 1;
      Options.plotColor  = gui_var.Types.plotColors(gui_var.runCount);
      Options.delta0     = RunSet(gui_var.runCount).RunStats.delta;
      RunStats           = RunSet(gui_var.runCount).RunStats;
      set(gui.runStatus,'String',['Resuming Run # ',int2str(gui_var.runCount)]);
      [BestF,BestI,RunStats,RunSet(1).Cache] = mads(Problem,iterate0, ...
                                                    Options,RunStats);
      delete(fullfile(gui_var.path,Problem.File.C));
      set([gui.CacheMenu.load; gui.CacheMenu.count; gui.CacheMenu.delete], ...
          'Enable','off');
   else
      set(gui.runStatus,'String', ...
                       ['Processing Run # ',int2str(gui_var.runCount+1)]);
      [BestF,BestI,RunStats,RunSet(1).Cache] = mads(Problem,iterate0,Options);
      gui_var.runCount = gui_var.runCount+1;
   end
   RunSet(gui_var.runCount).BestF    = BestF;
   RunSet(gui_var.runCount).BestI    = BestI;
   RunSet(gui_var.runCount).RunStats = RunStats;
   if (RunStats.time > 60)
      load train;
      sound(y);
   end

% Perform these tasks if error in MADS run
try
catch
   set(gui.runStatus,'String',['Run # ',int2str(gui_var.runCount+1),' failed']);
   set([gui.axesFilter; get(gui.axesFilter, 'Children')],'Visible','off');
   if (gui_var.runCount == 0)
      set([gui.axesHistory; get(gui.axesHistory,'Children')],'Visible','off');
   else
      set(gui.axesHistory,'ButtonDownFcn',gui_func.copyPlot);
   end
   RunSet(1).Cache = [];
   errordlg(lasterr,'MADS Runtime Error','modal'); beep
   rethrow(lasterror)
   rmappdata(0,'CACHE');
end

% Change figure window, as appropriate
set(gui.axesHistory, 'NextPlot','add');
set(gui.axesHistory, 'ButtonDownFcn',gui_func.copyPlot);
set(gui.axesFilter,  'ButtonDownFcn',gui_func.copyPlot);

if (gui_var.runCount > 0)
   set([gui.CacheMenu.save;gui.RunMenu.resume;gui.RunMenu.restart],'Enable','on');
   set([gui.Menu.Results; gui.ResultsMenu(gui_var.runCount)],'Visible','on');
   set(gui.runStatus,'String',['Run # ',int2str(gui_var.runCount),' complete']);
end
set(gui.stopRun,   'Visible','off','UserData',0);
set(gui.resumeRun, 'Visible','on');
setptr(fig,'arrow');
setappdata(fig,'gui_var',gui_var);
setappdata(fig,'RunSet',RunSet);

return

%*******************************************************************************
% loadMADS:  Assign GUI input fields to MADS variables.
% ------------------------------------------------------------------------------
% Called by: saveSession, runMADS
% Calls:     nomadm_compile
% VARIABLES:
%  Problem              = structure containing optimization problem data
%    .File              =   structure of problem file names
%      .F               =   name of functions file
%      .O               =   name of Omega file
%      .I               =   name of initial points file
%      .N               =   name of discrete neighbors file
%      .C               =   name of Cache File
%      .nameCache       =   name of the base workspace Cache variable
%      .fType           =   type of functions file (M=MATLAB,F=FORTRAN,C=C)
%  Options              = structure containing MADS parameters
%    .nSearches         =   number of Search types used
%    .Search(n)         =   structure of Search parameters
%      .type            =     string identifying the type of Search
%    .pollStrategy      =   string identifying selected Poll strategy
%    .pollOrder         =   string identifying selected Poll order strategy
%    .pollCenter        =   integer identifying selected Poll center
%    .pollComplete      =   turns on/off complete Polling
%    .loadCache         =   flag for loading a pre-existing Cache file
%    .countCache        =   flag for counting Cache points as function calls
%    .useFilter         =   use filter for nonlinear constraints
%    .removeRedundancy  =   remove redundant linear constraints
%    .runStochastic     =   flag for running as a stochastic problem
%    .accelerate        =   flag for accelerating mesh refinement
%    .scale             =   flag for scaling mesh directions
%    .plotHistory1      =   turns on/off a history plot
%    .plotHistory2      =   turns on/off a real-time history plot
%    .plotFilter        =   turns on/off a real-time filter plot
%    .runOneIteration   =   flag for running one MADS iteration at a time
%    .runUntilFeasible  =   flag for running MADS only until feasible
%    .runCount          =   MADS run counter
%    .hplothandle       =   handle for history plot axes
%    .fplothandle       =   handle for filter plot axes
%    .stophandle        =   handle for Stop Run pushbutton
%    .delta0            =   initial mesh size
%    .deltaMax          =   maximum mesh size
%    .meshRefine        =   mesh refinement factor
%    .meshCoarsen       =   mesh coarsening factor
%    .tolCache          =   tolerance for flagging point as being in Cache
%    .hmin              =   minimum h-value of an infeasible point
%    .hmax              =   maximum h-value of a filter point
%    .ePollTriggerF     =   f-value Extended Poll trigger
%    .ePollTriggerH     =   h-value Extended Poll trigger
%    .Term              =   substructure containing MADS termination criteria
%      .delta           =     mesh size parameter
%      .iter            =     maximum number of MADS iterations
%      .func            =     maximum number of function evaluations
%      .time            =     maximum CPU time
%      .fails           =     maximum number of consecutive Poll failures
%    .TermFlag          =   substructure of termination criteria on/off switches
%      .iter            =     turns on/off number of iterations
%      .nFunc           =     turns on/off number of function evaluations
%      .time            =     turns on/off CPU time
%      .nFails          =     turns on/off number of consecutive Poll failures
%      .relative        =     computes termination delta relative to .delta0
%  gui_var              = structure containing all GUI variables
%    .path              =   path of current optimization problem
%    .problemExt        =   filename extension of optimization problem 
%    .Options           =   current Options values
%    .Types             = lists of possible type
%      .Search          =   list of possible Search types
%      .poll            =   list of possible Poll strategies
%      .pollOrder       =   list of possible Poll order strategies
%    .Choice            = user choices
%      .pollStrategy    =   selected Poll strategy
%      .pollCenter      =   selected Poll center
%      .pollOrder       =   selected Poll order strategy
%  gui                  = structure of GUI object handles
%    .runStatus         =   handle for GUI figure window status bar
%    .problem           =   name of current optimization problem
%    .CacheMenu         =   handles for Cache menu items
%    .MADSMenu          =   handles for MADS menu items
%    .OptionsMenu       =   handles for Options menu items
%    .RunMenu           =   handles for Run menu items
%    .delta0            =   current initial mesh size
%    .deltaMax          =   current maximum mesh size
%    .meshRefine        =   current mesh refinement factor
%    .meshCoarsen       =   current mesh coarsening factor
%    .tolCache          =   current Cache tolerance
%    .hmin              =   current minimum infeasible h-value
%    .hmax              =   current maximum filter h-value
%    .ePollXiF          =   current f-value Extended Poll trigger
%    .ePollXiH          =   current h-value Extended Poll trigger
%    .Term              =   structure of handles for current termination criteria
%      .delta           =     current mesh size
%      .nIter           =     current maximum number of MADS iterations
%      .nFunc           =     current maximum number of function evaluations
%      .time            =     current maximum CPU time
%      .nFails          =   current max number of consecutive Poll failures
%    .TermFlag          =   structure of handles for termination checkboxes
%      .nIter           =     handle for number of iterations
%      .nFunc           =     handle for number of function evaluations
%      .time            =     handle for CPU time
%      .nFails          =     handle for number of consecutive Poll failures
%  nameProblem          = name of the optimization problem to be solved
%  language             = programming language of functions file
%  k                    = Search counter
%  field                = cell array of field names of the gui.Term substructure
%*******************************************************************************
function [Problem,Options] = loadMADS(gui,gui_var)

% Transfer Optimization Problem data from GUI into MADS input variables
addpath(gui_var.path);
Problem.nameCache = gui_var.nameCache;
nameProblem = get(gui.problem,'String');

ext = fieldnames(gui_var.FileExt);
for k = 1:length(ext)
   Problem.File.(ext{k}) = [nameProblem,gui_var.FileExt.(ext{k})];
end

% Compile non-Matlab Functions file, if possible
Problem.fType  = upper(gui_var.problemExt(2));
if ~strcmp(Problem.fType,'M')
   language = nomadm_compile(Problem.fType, ...
                           gui_var.path,Problem.File.F,gui_var.problemExt);
   set(gui.runStatus,'String',['Compiling ',language,' function file']);
end

% Transfer user options from GUI into MADS input variables
umtoggle(gui.scaleMenu(2));
umtoggle(gui.scaleMenu(3));
umtoggle(gui.CacheMenu.load);
umtoggle(gui.CacheMenu.count);
umtoggle(gui.MADSMenu.pollComplete);
umtoggle(gui.OptionsMenu.useFilter);
umtoggle(gui.OptionsMenu.removeRedundancy);
umtoggle(gui.OptionsMenu.runStochastic);
umtoggle(gui.OptionsMenu.accelerate);
umtoggle(gui.OptionsMenu.plotFilter);
umtoggle(gui.OptionsMenu.plotHistory1);
umtoggle(gui.OptionsMenu.plotHistory2);
umtoggle(gui.RunMenu.execFeasible);
umtoggle(gui.RunMenu.oneIteration);
umtoggle(gui.OptionsMenu.TermFlag.relative);

Options = gui_var.Options;
for k = 1:Options.nSearches
   Options.Search(k).type = gui_var.Types.search{gui_var.Choice.search(k)};
end
Options.pollStrategy      = gui_var.Types.poll{gui_var.Choice.pollStrategy};
Options.pollOrder         = gui_var.Types.pollOrder{gui_var.Choice.pollOrder};
Options.pollCenter        = gui_var.Choice.pollCenter - 1;
Options.pollComplete      = umtoggle(gui.MADSMenu.pollComplete);
Options.loadCache         = umtoggle(gui.CacheMenu.load);
Options.countCache        = umtoggle(gui.CacheMenu.count);
Options.useFilter         = umtoggle(gui.OptionsMenu.useFilter);
Options.removeRedundancy  = umtoggle(gui.OptionsMenu.removeRedundancy);
Options.runStochastic     = umtoggle(gui.OptionsMenu.runStochastic);
Options.accelerate        = umtoggle(gui.OptionsMenu.accelerate);
Options.scale             =  2*umtoggle(gui.scaleMenu(2)) + ...
                            10*umtoggle(gui.scaleMenu(3));
Options.plotFilter        = umtoggle(gui.OptionsMenu.plotFilter);
Options.plotHistory1      = umtoggle(gui.OptionsMenu.plotHistory1);
Options.plotHistory2      = umtoggle(gui.OptionsMenu.plotHistory2);
Options.runOneIteration   = umtoggle(gui.RunMenu.oneIteration);
Options.runUntilFeasible  = umtoggle(gui.RunMenu.execFeasible);
Options.TermFlag.relative = umtoggle(gui.OptionsMenu.TermFlag.relative);
Options.plotColor         = gui_var.Types.plotColors(gui_var.runCount+1);
Options.hplothandle       = gui.axesHistory;
Options.fplothandle       = gui.axesFilter;
Options.stophandle        = gui.stopRun;
Options.delta0            = str2double(get(gui.delta0,      'String'));
Options.deltaMax          = str2double(get(gui.deltaMax,    'String'));
Options.meshRefine        = str2double(get(gui.meshRefine,  'String'));
Options.meshCoarsen       = str2double(get(gui.meshCoarsen, 'String'));
Options.tolCache          = str2double(get(gui.tolCache,    'String'));
Options.hmin              = str2double(get(gui.hmin,        'String'));
Options.hmax              = str2double(get(gui.hmax,        'String'));
Options.ePollTriggerF     = str2double(get(gui.ePollXiF,    'String'));
Options.ePollTriggerH     = str2double(get(gui.ePollXiH,    'String'));

% Process Termination Criteria
field = fieldnames(gui.Term);
for k = 1:length(field)
   Options.Term.(field{k})     = str2double(get(gui.Term.(field{k}),'String'));
   Options.TermFlag.(field{k}) = get(gui.TermFlag.(field{k}),       'Value');
   if ~Options.TermFlag.(field{k}), Options.Term.(field{k}) = Inf; end
end
if (Options.TermFlag.relative)
   Options.Term.delta = Options.Term.delta*Options.delta0;
end

% Process RS Parameters
field = fieldnames(gui.RS);
for k = 1:length(field)
   Options.RS.(field{k}) = str2double(get(gui.RS.(field{k}),'String'));
end

return

%*******************************************************************************
% search_gui:  Displays a user input screen to set Search parameters.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: MADS-->Select Search Strategies)
% Calls:     modifySearchScreen, loadUserSearchFile, loadSearchOptions
% VARIABLES:
%  guiSearch         = structure for the Search GUI and its object handles
%    .fig            =   handle for Search figure window
%    .n              =   handle for Number of Search Types popup menu
%    .SurOptimizer   =   handle for surrogate optimizer popup menu
%    .Screen         = structure of handles for each k of n Searches
%      .label        =   handle for Search text label
%      .type         =   handle for popupmenu of possible Search types
%      .nIter        =   handle for number of iterations field
%      .nPoints      =   handle for number of Search points field
%      .file         =   handle for string containing name of Search file
%    .done           =   handle for Done pushbutton
%    .cancel         =   handle for Cancel pushbutton   
%  gui_var           = structure of GUI variables
%    .maxSearches    = maximum number of Search types that can be selected
%    .Labels         = long text labels used in popup menus
%      .Search       =   labels used in each Search type popup menu
%      .optimizer    =   labels for the surrogate optimizer to be used
%    .Choice.Search  = integers recording Search type popup menu choices
%    .Options.Search = vector of structures for each k of n Searches
%      .nIter        =   number of iterations
%      .nPoints      =   number of Search points
%      .file         =   string containing name of optional Search file
%  k                 = Search type counter
%  row               = row location of current Search figure window object
%  pathStr           = temporary storage for Search file path
%  filename          = temporary storage for Search filename
%  fileExt           = temporary storage for Search filename extension
%*******************************************************************************
function search_gui(h,event)

gui_var  = getappdata(gcbf,'gui_var');
gui_func = getappdata(gcbf,'gui_func');

%Set up "Search Options" figure window
guiSearch.fig = figure(...
   'Name',                           'Set Options for MADS SEARCH step', ...
   'Tag',                            'Search_GUI', ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.375, ...,
   'DefaultUIControlStyle',           'text', ...
   'WindowStyle',                     'modal', ...
   'Units',                           'normalized', ...
   'Position',                        [0.15 0.1 0.8 0.7], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

% Set up panel for Top deisplays
panelColor = get(guiSearch.fig,'DefaultUIControlBackgroundColor');

if gui_var.newMatlab
   guiSearch.panel = uipanel('Parent',guiSearch.fig, ...
                             'Position',[0, .9, 1, .1], ...
                             'BorderType','beveledin', ...
                             'BackgroundColor',panelColor);

   % Display Number of Search Types and Surrogate Optimizer
   uicontrol(guiSearch.panel, ...
       'Style',           'text', ...
       'String',          'Number of Search Types: ', ...
       'Position',        [.12 0 .18 .6]);
   guiSearch.n = uicontrol(guiSearch.panel, ...
       'Style',           'popupmenu', ...
       'String',          {'0','1','2','3','4','5','6','7','8'}, ...
       'BackgroundColor', 'white', ...
       'Value',           gui_var.Options.nSearches+1, ...
       'Position',        [.31 .03 .05 .65], ...
       'Callback',        gui_func.modifySearchScreen);
   uicontrol(guiSearch.panel, ...
       'Style',           'text', ...
       'String',          'Surrogate Optimizer: ', ...
       'Position',        [.45 0 .15 .6]);
   guiSearch.SurOptimizer = uicontrol(guiSearch.panel, ...
       'Style',           'popupmenu', ...
       'String',          gui_var.Labels.optimizer, ...
       'BackgroundColor', 'white', ...
       'Value',           gui_var.Choice.optimizer, ...
       'Position',        [.61 .03 .33 .65]);

else
   
   % Additional figure window partitioning lines and panel
   uicontrol(guiSearch.fig, 'Style','frame','Position', [0, .896, 1, .004],...
       'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
   uicontrol(guiSearch.fig, 'Style','frame','Position', [0, .9,   1, .004],...
       'ForegroundColor','white','BackgroundColor','white');
   uicontrol(guiSearch.fig, 'Style','frame','Position',[0,  .9,  1, .1], ...
      'ForegroundColor',panelColor,'BackgroundColor',panelColor);

   % Display Number of Search Types and Surrogate Optimizer
   uicontrol(guiSearch.fig, 'Style', 'text', ...
       'String', 'Number of Search Types: ', ...
       'Position', [.13 .905 .20 .06]);
   guiSearch.n = uicontrol(guiSearch.fig, ...
       'Style',           'popupmenu', ...
       'String',          {'0','1','2','3','4','5','6','7','8'}, ...
       'BackgroundColor', 'white', ...
       'Value',           gui_var.Options.nSearches+1, ...
       'Position',        [.35 .9125 .05 .06], ...
       'Callback',        gui_func.modifySearchScreen);
   uicontrol(guiSearch.fig, 'Style', 'text', ...
       'String', 'Surrogate Optimizer: ', ...
       'Position', [.48 .905 .15 .06]);
   guiSearch.SurOptimizer = uicontrol(guiSearch.fig, ...
       'Style',           'popupmenu', ...
       'String',          gui_var.Labels.optimizer, ...
       'BackgroundColor', 'white', ...
       'Value',           gui_var.Choice.optimizer, ...
       'Position',        [.63 .9125 .32 .06]);
end


% Text headers for the Search parameters
set(guiSearch.fig,'DefaultUIControlBackgroundColor',[.8,.8,.8]);
uicontrol(guiSearch.fig, 'Style', 'text', ...
    'String', 'Search Strategy', ...
    'Position', [.12 .80 .31 .06]);
uicontrol(guiSearch.fig, 'Style', 'text', ...
    'String', 'Number of Iterations', ...
    'Position', [.45 .80 .15 .06]);
uicontrol(guiSearch.fig, 'Style', 'text', ...
    'String', 'Number of Points', ...
    'Position', [.62 .80 .16 .06]);
uicontrol(guiSearch.fig, 'Style', 'text', ...
    'String', 'User File', ...
    'Position', [.80 .80 .16 .06]);

% Main loop for each of n Searches
for k = 1:gui_var.maxSearches
   row = .76 - .08*(k-1);

   guiSearch.Screen(k).label = uicontrol(guiSearch.fig, ...
      'Style',           'text', ...
      'String',          ['Search #', int2str(k), ':'], ...
      'Position',        [.03, row-.0075, .08, .06]);

   % The data fields for each Search
   guiSearch.Screen(k).type = uicontrol(guiSearch.fig, ...
      'Style',           'popupmenu', ...
      'String',          gui_var.Labels.search, ...
      'BackgroundColor', 'white', ...
      'Value',           gui_var.Choice.search(k), ...
      'Position',        [.12, row, .31, .06], ...
      'UserData',        {k,gui_var.Choice.search(k), ...
                         gui_var.Options.Search(k).local, ...
                         gui_var.Options.Search(k).merit}, ...
      'Callback',        gui_func.loadUserSearchFile);
   guiSearch.Screen(k).nIter = uicontrol(guiSearch.fig, ...
      'Style',           'edit', ...
      'String',          int2str(gui_var.Options.Search(k).nIter), ...
      'BackgroundColor', 'white', ...
      'FontSize',        .64, ...
      'Position',        [.45, row+.02, .15, .04]);
   guiSearch.Screen(k).nPoints = uicontrol(guiSearch.fig, ...
      'Style',           'edit', ...
      'String',          int2str(gui_var.Options.Search(k).nPoints), ...
      'BackgroundColor', 'white', ...
      'FontSize',        .64, ...
      'Position',        [.62, row+.02, .16, .04]);
   [pathStr,filename,fileExt] = fileparts(gui_var.Options.Search(k).file);
   guiSearch.Screen(k).file = uicontrol(guiSearch.fig, ...
      'Style',           'text', ...
      'String',          [filename,fileExt], ...
      'ForegroundColor', 'red', ...
      'Position',        [.80, row, .16, .06], ...
      'UserData',        pathStr);
end

% The Done and Cancel Buttons
guiSearch.done = uicontrol(guiSearch.fig, ...
   'Style',      'pushbutton', ...
   'String',     'Done', ...
   'FontWeight', 'bold', ...
   'Position',   [.33, .04, .10, .08], ...
   'UserData',   1, ...
   'Callback',   gui_func.loadSearchOptions);
guiSearch.cancel = uicontrol(guiSearch.fig, ...
   'Style',      'pushbutton', ...
   'String',     'Cancel', ...
   'FontWeight', 'bold', ...
   'Position',   [.50, .04, .10, .08], ...
   'UserData',   0, ...
   'Callback',   gui_func.loadSearchOptions);

setappdata(guiSearch.fig,'guiSearch',guiSearch);
modifySearchScreen(2,0);
return

%*******************************************************************************
% modifySearchScreen: Modify the Search Input Screen.
% ------------------------------------------------------------------------------
% Called by: search_gui (callback: Number of Search Types popupmenu
%            search_gui (directly)
% VARIABLES:
%  fig            = temporary storage of figure handles
%  gui_var        = structure of all GUI variables
%    .nSearches   =   handle for Number of Search Types popup menu
%    .maxSearches =   maximum number of Search types allowed
%  guiSearch      = structure of all Search window object handles
%    .n           =   number of Search types popup menu
%    .Screen(k)   =   Search fields for the Search #k
%      .label     =     text label appearing on the Search Screen
%      .type      =     type of Search
%      .nIter     =     number of iterations
%      .nPoints   =     number of Search points
%      .file      =     file used during Search
%    .done        =   handle for Search input screen Done button
%    .cancel      =   handle for Search Input Screen Cancel button   
%  nSearches      = number of Search types selected
%  pos            = screen position of the final Search field
%*******************************************************************************
function modifySearchScreen(h,event)

% Get application data
fig       = findobj('Tag','NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');
fig       = findobj('Tag','Search_GUI');
guiSearch = getappdata(fig,'guiSearch');

% Turn on/off appropriate Search fields
nSearches = get(guiSearch.n,'Value') - 1;
for k = 1:gui_var.maxSearches
   if (k <= nSearches)
      set(guiSearch.Screen(k).label,   'Visible','on');
      set(guiSearch.Screen(k).type,    'Visible','on');
      set(guiSearch.Screen(k).nIter,   'Visible','on');
      set(guiSearch.Screen(k).nPoints, 'Visible','on');
      set(guiSearch.Screen(k).file,    'Visible','on');
   else
      set(guiSearch.Screen(k).label,   'Visible','off');
      set(guiSearch.Screen(k).type,    'Visible','off','Value', 1);
      set(guiSearch.Screen(k).nIter,   'Visible','off','String','1');
      set(guiSearch.Screen(k).nPoints, 'Visible','off','String','1');
      set(guiSearch.Screen(k).file,    'Visible','off','String','');
   end
end

% Move Done and Cancel buttons
if nSearches
   pos = get(guiSearch.Screen(nSearches).type,'Position');
else
   pos = [.01, .82];
end
set(guiSearch.done,   'Position', [.33, pos(2)-.16, .10, .08]);
set(guiSearch.cancel, 'Position', [.50, pos(2)-.16, .10, .08]);
return

%*******************************************************************************
% loadUserSearchFile: Load a user-specified Search file.
% ------------------------------------------------------------------------------
% Called by: search_gui (callback: each Search Strategy popup menu)
% Calls:     dace_gui, nw_gui
% VARIABLES:
%  fig                   = temporary storage of figure handles
%  gui_var               = structure of all GUI variables
%    .Types.search       =   list of possible Search types
%    .maxSearches        =   maximum number of Searches that can be done
%  guiSearch             = handles for all Search window objects
%    .Screen.file        =   handle for Search figure window user file
%  newValue              = number associated with Search menu choice
%  userData              = temporary storage of GUI object user data
%  k                     = Search number
%  previousType          = previously selected Search type
%  gui_func              = handles of all NOMADm callback functions
%    .loadUserSearchFile = handle for the function of the same name
%  SName, SPath          = name and path of user Search file
%  spec                  = file types for use in input dialog boxes
%*******************************************************************************
function loadUserSearchFile(h,event)

% Get application data
fig       = findobj('Tag','NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');
gui_func  = getappdata(fig,'gui_func');
fig       = findobj('Tag','Search_GUI');
guiSearch = getappdata(fig,'guiSearch');

newValue  = get(gcbo, 'Value');
userData  = get(gcbo, 'UserData');
[k,previousType,local,merit] = deal(userData{:});

switch gui_var.Types.search{newValue}

% Load a custom Search or Surrogate file 
case {'Custom','CustomS'}
   spec = {'*.m',        'Matlab M-files (*.m)'; ...
           '*.f; *.for', 'Fortran files (*.f,*.for)'; ...
           '*.c; *.C',   'C/C++ files (*.c,*.C)'};
   [sName,sPath] = uigetfile(spec,['Choose File for Search #',int2str(k)]);
   if (sName)
      set(guiSearch.Screen(k).file, 'UserData', sPath);
      set(guiSearch.Screen(k).file, 'String',   sName);
      set(gcbo, 'UserData',{k,newValue});
   else
      set(gcbo, 'Value', previousType);
   end

% Load a DACE surrogate
case {'DACE'}
   dace_gui(k,gui_var,gui_func);
   
% Load a NW surrogate
case {'NW'}
   nw_gui(k,gui_var,gui_func);

% Do not load a file
otherwise
   set(guiSearch.Screen(k).file, 'String', '');
   set(gcbo, 'UserData',{k,newValue,local,merit});
end
return

%*******************************************************************************
% loadSearchOptions: Load selected Search parameters.
% ------------------------------------------------------------------------------
% Called by: search_gui (callback: Done and Cancel pushbuttons
% Calls:     updateSearchLabels
% VARIABLES:
%  k                = Search number
%  fig              = temporary storage of the appropriate figure handle
%  gui_var          = structure of all GUI handles and variables
%    .Choice        =   structure of user choices
%      .search(k)   =     user Search choices
%      .optimizer   =     choice of surrogate optimizer
%    .Options       =   structure of MADS parameters settings
%      .Search(k)   =     structure of user-selected k-th Search
%      .nSearches   =     number of Search types to be used
%    .Types         =   lists of possible types
%      .Search      =   list of possible Search types
%      .optimizer   =   list of possible surrogate optimizers
%    .noGrad        =   flag indicating initial availability of gradients
%  guiSearch        = handles for all Search window objects
%    .n             = handle for the Number of Searches field
%    .SurOptimizer = handle for the Surrogate Optimizer popup menu
%    .Screen        = object handles for each Search
%      .type        =   string identifying Search type
%      .label       =   long text label for Search type
%      .nIter       =   number of iterations to perform Search
%      .nPoints     =   number of Search points
%      .file        =   optional user file defining Search
%  loadSearch       = flag indicating if Search options will be loaded
%  Search           = temporary storage of gui_var.Options.Search(k)
%  maxDisplay       = number of Search Types displayed on the main GUI
%  searchLabel      = Search label that appears on the main GUI
%*******************************************************************************
function loadSearchOptions(h,event)

fig       = findobj('Tag','NOMADm_GUI');
gui       = getappdata(fig,'gui');
gui_var   = getappdata(fig,'gui_var');
guiSearch = getappdata(gcbf,'guiSearch');

loadSearch = get(gcbo,'UserData');
if loadSearch
   gui_var.Options.nSearches    = get(guiSearch.n,'Value') - 1;
   gui_var.Choice.optimizer     = get(guiSearch.SurOptimizer,'Value');
   gui_var.Options.SurOptimizer = ...
                   gui_var.Types.optimizer{gui_var.Choice.optimizer};
   for k = 1:gui_var.Options.nSearches
      gui_var.Choice.search(k) = get(guiSearch.Screen(k).type, 'Value');
      Search.type = gui_var.Types.search{gui_var.Choice.search(k)};
      if (gui_var.noGrad && strcmp(Search.type, 'GPollI'))
         uiwait(msgbox('No derivatives available for this problem', ...
                      ['Error in Search Type #', int2str(k)],'error','modal'));
         return
      end
      Search.label   = gui_var.Labels.search{gui_var.Choice.search(k)};
      Search.nIter   = str2double(get(guiSearch.Screen(k).nIter,  'String'));
      Search.nPoints = str2double(get(guiSearch.Screen(k).nPoints,'String'));
      Search.file    = fullfile(...
                          get(guiSearch.Screen(k).file, 'UserData'), ...   
                          get(guiSearch.Screen(k).file, 'String'));
      UserData = get(guiSearch.Screen(k).type, 'UserData');
      Search.local   = UserData{3};
      Search.merit   = UserData{4};
      gui_var.Options.Search(k) = Search;
   end

   % Update GUI Search fields
   maxDisplay  = length(gui.searchLabel);
   searchLabel = updateSearchLabels(maxDisplay, ...
                       gui_var.Options.nSearches, gui_var.Options.Search);
   for k = 1:maxDisplay
      set(gui.searchLabel(k), 'String', searchLabel{k});
   end
end

setappdata(fig,'gui_var',gui_var);
close(guiSearch.fig);
return

%*******************************************************************************
% dace_gui:  Displays a user input screen to set DACE Toolbox parameters.
% ------------------------------------------------------------------------------
% Called by: loadUserSearchFile
% Calls:     loadDACEOptions
% VARIABLES:
%  k                    = Search number for this DACE screen
%  guiDACE              = structure of all DACE window object handles
%    .fig               =   DACE figure window
%    .daceRegression    =   popup menu of regression functions
%    .daceCorrelation   =   popup menu of correlation functions
%    .daceTheta         =   field for estimating correlation parameter
%    .daceLower         =   field for entering lower bound for theta
%    .daceUpper         =   field for entering upper bound for theta
%    .daceIsotropic     =   checkbox for isotropic correlations
%    .done              =   the Done pushbutton
%    .cancel            =   the Cancel pushbutton
%  fig                  = handle for the NOMADm figure window
%  gui_var              = structure of all GUI variables
%    .Labels            =   labels for DACE function popup menus
%      .daceRegression  =     labels for DACE regression functions
%      .daceCorrelation =     labels for DACE correlation functions
%    .Choice            =   integer choices for DACE functions
%      .daceReg         =     choice for DACE regression function
%      .daceCorr        =     choice for DACE correlation function
%  gui_func             = structure of all GUI function handles
%  Options              = structure of MADS options
%    .dace              =   user-chosen DACE Toolbox parameters
%      .theta           =   estimate for theta
%      .lower           =   lower bound for theta
%      .upper           =   upper bound for theta
%      .isotropic       =   flag for isotropic theta
%*******************************************************************************
function guiDACE = dace_gui(k,gui_var,gui_func)

% Set up "DACE Options" figure window
guiDACE.fig = figure(...
   'Name', ['DACE Toolbox Options (Search #', int2str(k),')'], ...
   'Tag',  'DACE_GUI', ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.5, ...,
   'DefaultUIControlStyle',           'text', ...
   'WindowStyle',                     'modal', ...
   'Resize',                          'on', ...
   'Units',                           'normalized', ...
   'Position',                        [0.25 0.30 0.40 0.40], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');
if gui_var.newMatlab
   uipanel;
else
   panelColor = get(guiDACE.fig,'DefaultUIControlBackgroundColor');
   uicontrol(guiDACE.fig, 'Style','frame','Position',[0, 0, 1, 1], ...
      'ForegroundColor',panelColor,'BackgroundColor',panelColor);

end

% Figure window shading to make it appear recessed
uicontrol(guiDACE.fig, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiDACE.fig, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiDACE.fig, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(guiDACE.fig, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for DACE screen objects and data fields
uicontrol(guiDACE.fig, ...
      'Style',               'text', ...
      'String',              'Regression Model: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .865 .34 .08]);
guiDACE.reg                = uicontrol(guiDACE.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.daceRegression, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.daceReg(k), ...
      'Position',            [.37 .87 .60 .09]);

uicontrol(guiDACE.fig, ...
      'Style',               'text', ...
      'String',              'Correlation Model: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .775 .34 .08]);
guiDACE.corr               = uicontrol(guiDACE.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.daceCorrelation, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.daceCorr(k), ...
      'Position',            [.37 .78 .60 .09]);

uicontrol(guiDACE.fig, ...
      'Style',               'text', ...
      'String',              'Theta Initial Guess: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .625 .34 .08]);
guiDACE.theta              = uicontrol(guiDACE.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).theta), ...
      'FontSize',            .6, ...
      'BackgroundColor',     'white', ...
      'Position',            [.37 .64 .35 .08]);

uicontrol(guiDACE.fig, ...
      'Style',               'text', ...
      'String',              'Theta Lower Bound: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .535 .34 .08]);
guiDACE.lower              = uicontrol(guiDACE.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).lower), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .55 .35 .08]);

uicontrol(guiDACE.fig, ...
      'Style',               'text', ...
      'String',              'Theta Upper Bound: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .445 .34 .08]);
guiDACE.upper              = uicontrol(guiDACE.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.dace(k).upper), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .46 .35 .08]);

guiDACE.isotropic          = uicontrol(guiDACE.fig, ...
      'Style',               'checkbox',  ...
      'String',              'Isotropic', ...
      'Value',               gui_var.Options.dace(k).isotropic,   ...
      'Position',            [.75 .55 .20 .09]);
  
guiDACE.local              = uicontrol(guiDACE.fig, ...
      'Style',               'checkbox', ...
      'String',              'Restrict Model to Trust Region', ...
      'Value',               gui_var.Options.Search(k).local, ...
      'Position',            [.02 .30 .70 .09]);
guiDACE.merit              = uicontrol(guiDACE.fig, ...
      'Style',               'checkbox', ...
      'String',              'Use Merit Function to Penalize Clustering', ...
      'Value',               gui_var.Options.Search(k).merit, ...
      'Position',            [.02 .21 .70 .09]);

% The Done and Cancel Buttons
guiDACE.done                 = uicontrol(guiDACE.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Done', ...
      'FontWeight',          'bold', ...
      'Position',            [.27, .02, .18, .12], ...
      'Callback',            {gui_func.loadDACEOptions, 1, k});
guiDACE.cancel             = uicontrol(guiDACE.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Cancel', ...
      'FontWeight',          'bold', ...
      'Position',            [.54, .02, .18, .12], ...
      'Callback',            {gui_func.loadDACEOptions, 0, k});
setappdata(guiDACE.fig,'guiDACE',guiDACE);
return

%*******************************************************************************
% loadDACEOptions: Select a DACE surrogate option.
% ------------------------------------------------------------------------------
% Called by: dace_gui (callback: Done and Cancel pushbuttons)
% VARIABLES:
%  loadDACE           = flag indicating if DACE options are accepted
%  k                  = Search number
%  figSearch          = handle for Search figure window
%  guiSearch          = structure of Search window object handles
%    .Screen          =   handles for Search figure window objects
%      .type          =     type of Search
%      .file          =     optional user file defining Search
%    .daceXXX         =   handles for DACE figure window objects
%  figDACE            = handle for DACE figure window
%  guiDACE            = structure of handles for DACE window objects
%    .Reg             =   handle for regression function popup menu
%    .Corr            =   handle for correlation function popup menu
%    .Theta           =   handle for correlation parameter initial guess
%    .Lower           =   handle for correlation parameter lower bound
%    .Upper           =   handle for correlation parameter upper bound
%    .Isotropic       =   handle for isotropic parameters check box
%  fig                = handle for the NOMADm figure window
%  gui_var            = structure of all GUI variables
%    .Types           =   lists of possible types
%      .daceReg       =     list of possible DACE regression functions
%      .daceCorr      =     list of possible DACE correlation functions
%    .Options.dace(k) =   structure of DACE parameters for Search k
%  previousType       = previously selected Search type
%  dace               = structure of DACE parameters
%    .Reg             =   regression function handle
%    .Corr            =   correlation function handle
%    .Theta           =   initial guess for correlation parameters
%    .Lower           =   lower bounds for correlation parameters
%    .Upper           =   upper bounds for correlation parameters
%    .Isotropic       =   flag for isotropic correlation parameters
%*******************************************************************************
function loadDACEOptions(h,event,loadDACE,k)

figSearch = findobj('Tag','Search_GUI');
guiSearch = getappdata(figSearch,'guiSearch');
figDACE   = findobj('Tag','DACE_GUI');
guiDACE   = getappdata(figDACE,'guiDACE');
fig       = findobj('Tag', 'NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');

previousType = get(guiSearch.Screen(k).type,'UserData');
if loadDACE
   dace.reg       = gui_var.Types.daceReg{ get(guiDACE.reg, 'Value')};
   dace.corr      = gui_var.Types.daceCorr{get(guiDACE.corr,'Value')};
   dace.theta     = str2double(get(guiDACE.theta,'String'));
   dace.lower     = str2double(get(guiDACE.lower,'String'));
   dace.upper     = str2double(get(guiDACE.upper,'String'));
   dace.isotropic = get(guiDACE.isotropic,'Value');
   gui_var.Options.dace(k) = dace;
   set(guiSearch.Screen(k).file, 'UserData', '','String',dace.reg);
   set(guiSearch.Screen(k).type, ...
      'UserData', {k,get(guiSearch.Screen(k).type,'Value'), ...
                     get(guiDACE.local,'Value'), get(guiDACE.merit,'Value')});
else
   set(guiSearch.Screen(k).type, 'Value', previousType{2});
end
close(guiDACE.fig);
return

%*******************************************************************************
% nw_gui:  Displays a user input screen to set NW Toolbox parameters.
% ------------------------------------------------------------------------------
% Called by: loadUserSearchFile
% Calls:     loadNWOptions
% VARIABLES:
%  k                    = Search number for this NW screen
%  guiNW                = structure of all NW window object handles
%    .fig               =   NW figure window
%    .nwKernel          =   popup menu of kernel functions
%    .nwLocal           =   checkbox for restricting model to local region
%    .done              =   the Done pushbutton
%    .cancel            =   the Cancel pushbutton
%  fig                  = handle for the NOMADm figure window
%  gui_var              = structure of all GUI variables
%    .Labels            =   labels for NW function popup menus
%      .nwKernel        =     labels for NW Kernel functions
%    .Choice            =   integer choices for NW functions
%      .nwKernel        =     choice for NW kernel function
%  gui_func             = structure of all GUI function handles
%  Options              = structure of MADS options
%    .nw                =   user-chosen NW Toolbox parameters
%      .kernel          =   estimate for theta
%*******************************************************************************
function guiNW = nw_gui(k,gui_var,gui_func)

% Set up "NW Options" figure window
guiNW.fig = figure(...
   'Name', ['NW Toolbox Options (Search #', int2str(k),')'], ...
   'Tag',  'NW_GUI', ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.5, ...,
   'DefaultUIControlStyle',           'text', ...
   'WindowStyle',                     'modal', ...
   'Resize',                          'on', ...
   'Units',                           'normalized', ...
   'Position',                        [0.25 0.30 0.40 0.40], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');
if gui_var.newMatlab
   uipanel;
else
   panelColor = get(guiNW.fig,'DefaultUIControlBackgroundColor');
   uicontrol(guiNW.fig, 'Style','frame','Position',[0, 0, 1, 1], ...
      'ForegroundColor',panelColor,'BackgroundColor',panelColor);

end

% Figure window shading to make it appear recessed
uicontrol(guiNW.fig, 'Style','frame','Position', [0, 0, 1, .004],    ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiNW.fig, 'Style','frame','Position', [.997, 0, .003, 1], ...
          'ForegroundColor','white','BackgroundColor','white');
uicontrol(guiNW.fig, 'Style','frame','Position', [0, .996, 1, .004], ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
uicontrol(guiNW.fig, 'Style','frame','Position', [0, 0, .003, 1],    ...
          'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);

% Labels for NW screen objects and data fields
uicontrol(guiNW.fig, ...
      'Style',               'text', ...
      'String',              'Kernel Function: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .855 .34 .08]);
guiNW.kernel                = uicontrol(guiNW.fig, ...
      'Style',               'popupmenu', ...
      'String',              gui_var.Labels.nwKernel, ...
      'BackgroundColor',     'white', ...
      'Value',               gui_var.Choice.nwKernel(k), ...
      'Position',            [.37 .86 .60 .09]);

uicontrol(guiNW.fig, ...
      'Style',               'text', ...
      'String',              'Sigma Initial Guess: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .625 .34 .08]);
guiNW.sigma                = uicontrol(guiNW.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.nw(k).sigma), ...
      'FontSize',            .6, ...
      'BackgroundColor',     'white', ...
      'Position',            [.37 .64 .35 .08]);

uicontrol(guiNW.fig, ...
      'Style',               'text', ...
      'String',              'Sigma Lower Bound: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .535 .34 .08]);
guiNW.lower              = uicontrol(guiNW.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.nw(k).lower), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .55 .35 .08]);

uicontrol(guiNW.fig, ...
      'Style',               'text', ...
      'String',              'Sigma Upper Bound: ', ...
      'HorizontalAlignment', 'right', ...
      'FontSize',            .6, ...
      'Position',            [.02 .445 .34 .08]);
guiNW.upper              = uicontrol(guiNW.fig, ...
      'Style',               'edit', ...
      'String',              num2str(gui_var.Options.nw(k).upper), ...
      'BackgroundColor',     'white', ...
      'FontSize',            .6, ...
      'Position',            [.37 .46 .35 .08]);

guiNW.local                = uicontrol(guiNW.fig, ...
      'Style',               'checkbox', ...
      'String',              'Restrict Model to Trust Region', ...
      'Value',               gui_var.Options.Search(k).local, ...
      'Position',            [.02 .30 .70 .09]);
guiNW.merit                = uicontrol(guiNW.fig, ...
      'Style',               'checkbox', ...
      'String',              'Use Merit Function to Penalize Clustering', ...
      'Value',               gui_var.Options.Search(k).merit, ...
      'Position',            [.02 .21 .70 .09]);

% The Done and Cancel Buttons
guiNW.done                 = uicontrol(guiNW.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Done', ...
      'FontWeight',          'bold', ...
      'Position',            [.27, .02, .18, .12], ...
      'Callback',            {gui_func.loadNWOptions, 1, k});
guiNW.cancel             = uicontrol(guiNW.fig, ...
      'Style',               'pushbutton', ...
      'String',              'Cancel', ...
      'FontWeight',          'bold', ...
      'Position',            [.54, .02, .18, .12], ...
      'Callback',            {gui_func.loadNWOptions, 0, k});
setappdata(guiNW.fig,'guiNW',guiNW);
return

%*******************************************************************************
% loadNWOptions: Select a NW surrogate option.
% ------------------------------------------------------------------------------
% Called by: nw_gui (callback: Done and Cancel pushbuttons)
% VARIABLES:
%  loadNW             = flag indicating if NW options are accepted
%  k                  = Search number
%  figSearch          = handle for Search figure window
%  guiSearch          = structure of Search window object handles
%    .Screen          =   handles for Search figure window objects
%      .type          =     type of Search
%      .file          =     optional user file defining Search
%    .nwXXX           =   handles for NW figure window objects
%  figNW              = handle for NW figure window
%  guiNW              = structure of handles for NW window objects
%    .kernel          =   handle for kernel function popup menu
%    .local           =   checkbox handle for using a trust region
%  fig                = handle for the NOMADm figure window
%  gui_var            = structure of all GUI variables
%    .Types           =   lists of possible types
%      .nwKernel      =     list of possible NW regression functions
%    .Options.nw(k)   =   structure of NW parameters for Search k
%  previousType       = previously selected Search type
%  nw                 = structure of NW parameters
%    .kernel          =   kernel function handle
%*******************************************************************************
function loadNWOptions(h,event,loadNW,k)

figSearch = findobj('Tag','Search_GUI');
guiSearch = getappdata(figSearch,'guiSearch');
figNW     = findobj('Tag','NW_GUI');
guiNW     = getappdata(figNW,'guiNW');
fig       = findobj('Tag', 'NOMADm_GUI');
gui_var   = getappdata(fig,'gui_var');

previousType = get(guiSearch.Screen(k).type,'UserData');
if loadNW
   nw.kernel = gui_var.Types.nwKernel{get(guiNW.kernel,'Value')};
   nw.sigma  = str2double(get(guiNW.sigma,'String'));
   nw.lower  = str2double(get(guiNW.lower,'String'));;
   nw.upper  = str2double(get(guiNW.upper,'String'));;
   gui_var.Options.nw(k) = nw;
   
   set(guiSearch.Screen(k).file,'UserData','','String',nw.kernel);
   set(guiSearch.Screen(k).type, ...
      'UserData', {k,get(guiSearch.Screen(k).type,'Value'), ...
                   get(guiNW.local,'Value'), get(guiNW.merit,'Value')});
else
   set(guiSearch.Screen(k).type, 'Value', previousType{2});
end
close(guiNW.fig);
return

%*******************************************************************************
% updateScreen:  Displays a user input screen to set Search parameters.
% ------------------------------------------------------------------------------
% Called by: loadSession, resetSession
% Calls:     updateSearchLabels
% VARIABLES:
%  Choice               = structure of menu choices
%    .pollStrategy      =   integer choice of Poll strategy
%    .pollOrder         =   integer choice of Poll order strategy
%    .pollCenter        =   integer choice of Poll center
%  Options              = structure of MADS parameters
%    .nSearches         =   number of Search types used
%    .Search(n)         =   structure of Search parameters
%    .delta0            =   initial mesh size
%    .deltaMax          =   maximum mesh size
%    .meshRefine        =   mesh refinement factor
%    .meshCoarsen       =   mesh coarsening factor
%    .tolCache          =   tolerance for flagging point as being in Cache
%    .hmin              =   minimum h-value of an infeasible point
%    .hmax              =   maximum h-value of a filter point
%    .ePollTriggerF     =   f-value Extended Poll trigger
%    .ePollTriggerH     =   h-value Extended Poll trigger
%    .Term              =   substructure containing MADS termination criteria
%      .delta           =     mesh size parameter
%      .iter            =     maximum number of MADS iterations
%      .func            =     maximum number of function evaluations
%      .time            =     maximum CPU time
%      .fails           =     maximum number of consecutive Poll failures
%    .TermFlag          =   substructure of termination criteria on/off switches
%      .iter            =     turns on/off number of iterations
%      .nFunc           =     turns on/off number of function evaluations
%      .time            =     turns on/off CPU time
%      .nFails          =     turns on/off number of consecutive Poll failures
%    .scale             =   flag for scaling mesh directions
%    .pollComplete      =   turns on/off complete Polling
%    .TermFlag.relative =   computes termination delta relative to .delta0
%    .useFilter         =   use filter for nonlinear constraints
%    .removeRedundancy  =   remove redundant linear constraints
%    .runStochastic     =   run as stochastic optimization problem
%    .accelerate        =   flag for accelerating mesh refinement
%    .plotHistory1      =   turns on/off a history plot
%    .plotHistory2      =   turns on/off a real-time history plot
%    .plotFilter        =   turns on/off a real-time filter plot
%    .loadCache         =   flag for loading a pre-existing Cache file
%    .countCache        =   flag for counting Cache points as function calls
%    .runOneIteration   =   flag for running one MADS iteration at a time
%    .runUntilFeasible  =   flag for running MADS only until feasible
%  gui_var.Labels       = sub-structure of GUI labels
%    .pollStrategy      =   Poll strategy label
%    .pollCenter        =   Poll center label
%    .pollOrder         =   Poll order label
%  gui                  = structure of GUI object handles
%    .SearchType(k)     =   current Search types
%    .pollStrategy      =   current Poll strategy
%    .pollCenter        =   current Poll center
%    .pollOrder         =   current Poll order type
%    .delta0            =   current initial mesh size
%    .deltaMax          =   current maximum mesh size
%    .meshRefine        =   current mesh refinement factor
%    .meshCoarsen       =   current mesh coarsening factor
%    .tolCache          =   current Cache tolerance
%    .hmin              =   current minimum infeasible h-value
%    .hmax              =   current maximum filter h-value
%    .ePollXiF          =   current f-value Extended Poll trigger
%    .ePollXiH          =   current h-value Extended Poll trigger
%    .TermDelta         =   current mesh size termination criteria
%    .TermIter          =   current maximum number of MADS iterations
%    .TermFunc          =   current maximum number of function evaluations
%    .TermTime          =   current maximum CPU time
%    .TermFails         =   current max number of consec Poll failures
%    .TermFlag          =   structure of handles for termination checkboxes
%      .nIter           =     handle for number of iterations
%      .nFunc           =     handle for number of function evaluations
%      .time            =     handle for CPU time
%      .nFails          =     handle for number of consecutive Poll failures
%    .CacheMenu         =   Cache menu items
%    .MADSMenu          =   MADS menu items
%    .OptionsMenu       =   Options menu items
%    .RunMenu           =   Run menu items
%  onoff                = cell array of two strings "on" and "off"
%  maxDisplay           = number of Search Types displayed on the main GUI
%  searchLabel          = Search label that appears on the main GUI
%*******************************************************************************
function updateScreen(Choice,Options)

fig     = findobj('Tag','NOMADm_GUI');
gui_var = getappdata(fig,'gui_var');
gui     = getappdata(fig,'gui');
onoff   = {'on','off'};

% Update GUI Search fields
maxDisplay  = length(gui.searchLabel);
searchLabel = updateSearchLabels(maxDisplay,Options.nSearches,Options.Search);
for k = 1:maxDisplay
   set(gui.searchLabel(k), 'String', searchLabel{k});
end

% Update other GUI fields
set(gui.pollStrategy,'String',gui_var.Labels.pollStrategy(Choice.pollStrategy));
set(gui.pollCenter,  'String',gui_var.Labels.pollCenter(Choice.pollCenter));
set(gui.pollOrder,   'String',gui_var.Labels.pollOrder(Choice.pollOrder));
set(gui.delta0,      'String', num2str(Options.delta0,       '%1.1f'));
set(gui.deltaMax,    'String', num2str(Options.deltaMax,     '%2.1f'));
set(gui.meshRefine,  'String', num2str(Options.meshRefine,   '%1.1f'));
set(gui.meshCoarsen, 'String', num2str(Options.meshCoarsen,  '%1.1f'));
set(gui.tolCache,    'String', num2str(Options.tolCache,     '%1.6g'));
set(gui.hmin,        'String', num2str(Options.hmin,         '%1.5g'));
set(gui.hmax,        'String', num2str(Options.hmax,         '%1.2f'));
set(gui.ePollXiF,    'String', num2str(Options.ePollTriggerF,'%1.5g'));
set(gui.ePollXiH,    'String', num2str(Options.ePollTriggerH,'%1.5g'));
set(gui.Term.delta,  'String', num2str(Options.Term.delta,   '%1.5g'));
field = fieldnames(gui.TermFlag);
for k = 1:length(field)
   set(gui.Term.(field{k}),     'String', Options.Term.(field{k}));
   set(gui.TermFlag.(field{k}), 'Value',  Options.TermFlag.(field{k}));
end
field = fieldnames(gui.RS);
for k = 1:length(field)
   set(gui.RS.(field{k}), 'String', Options.RS.(field{k}));
end
set(gui.scaleMenu(1),    'Checked', onoff{1+~(Options.scale == 0)});
set(gui.scaleMenu(2),    'Checked', onoff{1+~(Options.scale == 2)});
set(gui.scaleMenu(3),    'Checked', onoff{1+~(Options.scale == 10)});
set(gui.CacheMenu.load,  'Checked', onoff{2-Options.loadCache});
set(gui.CacheMenu.count, 'Checked', onoff{2-Options.countCache});
set(gui.MADSMenu.pollComplete,        'Checked', ...
    onoff{2-Options.pollComplete});
set(gui.OptionsMenu.TermFlag.relative,'Checked', ...
    onoff{2-Options.TermFlag.relative});
set(gui.OptionsMenu.useFilter,        'Checked', ...
    onoff{2-Options.useFilter});
set(gui.OptionsMenu.removeRedundancy, 'Checked', ...
    onoff{2-Options.removeRedundancy});
set(gui.OptionsMenu.runStochastic,    'Checked', ...
    onoff{2-Options.runStochastic});
set(gui.OptionsMenu.accelerate,  'Checked', onoff{2-Options.accelerate});
set(gui.OptionsMenu.plotHistory1,'Checked', onoff{2-Options.plotHistory1});
set(gui.OptionsMenu.plotHistory2,'Checked', onoff{2-Options.plotHistory2});
set(gui.OptionsMenu.plotFilter,  'Checked', onoff{2-Options.plotFilter});
set(gui.RunMenu.oneIteration,    'Checked', onoff{2-Options.runOneIteration});
set(gui.RunMenu.execFeasible,    'Checked', onoff{2-Options.runUntilFeasible});
set(gui.MADSMenu.editRSParam,    'Enable',  onoff{2-Options.runStochastic});

gui_func = getappdata(gcf,'gui_func');
gui_func.toggleStochastic(1,0);
gui_func.toggleStochastic(1,0);
return

%*******************************************************************************
% updateSearchLabels:  Displays the Search Labels on the main GUI.
% ------------------------------------------------------------------------------
% Called by: loadSearchOptions, updateScreen
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
%*******************************************************************************
function label = updateSearchLabels(maxDisplay,nSearches,Search)

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
    case {'NW'}
       label{k} = ['N-W Surrogate: ',    searchFile];
    case {'Custom'}
       label{k} = ['Custom Search: ',    searchFile];
    case {'CustomS'}
       label{k} = ['Custom Surrogate: ', searchFile];
    end
    label{k} = [label{k}, ' (', int2str(Search(k).nIter), ')'];
end
return

%*******************************************************************************
% results_gui:  Displays Results from MADS run.
% ------------------------------------------------------------------------------
% Called by: nomadm_gui (callback: Results-->View Run #k)
% Calls:     solutions_gui
% VARIABLES:
%  k               = run number
%  gui_var         = structure of GUI variables
%  gui_func        = structure of NOMADm figure window function handles 
%  gui             = structure of NOMADm window object handles 
%    .problem      =   display of the optimization problem name
%    .Term         =   handles of termination criteria display
%      .delta      =     display of mesh size termination criteria
%      .nIter      =     display of maximum number of MADS iterations
%      .nFunc      =     display of maximum number of function evaluations
%      .time       =     display of maximum CPU time
%      .nFails     =     display of maximum number of consecutive Poll failures
%  RunSet          = structure containing all information about a MADS run
%    .BestF        =   best feasible iterate for Run k
%    .BestI        =   least infeasible iterate for Run k
%    .RunStats     =   MADS Run k statistics (iterations, CPU time, etc.)
%      .delta      =     final mesh size
%      .nIter      =     total number of iterations
%      .nFunc      =     total number of function evaluations
%      .time       =     total CPU time expended
%      .nFails     =     final number of consecutive Poll failures
%      .stopRun    =     final value of Stop Run button
%      .nCacheHits =     total number of Cache hits
%  k               = run number
%  digits          = number of digits accuracy to present on screen
%  noYes           = cell array containing the strings "No" and "Yes"
%  BestF, BestI    = best feasible and least infeasible solutions found
%  isRS            = flag indicating if the R&S method was applied
%  fig             = handle for the Results figure window
%  Run             = structure of handles for run statistics display
%    .delta        =   final mesh size
%    .nIter        =   number of iterations
%    .nFunc        =   number of function evaluations
%    .nGrad        =   number of gradient evaluations
%    .time         =   CPU time expended
%    .nFails       =   number of consecutive Poll failures
%    .stop         =   interrupted by user (Yes/No)
%    .nCacheHits   =   number of Cache hits
%    .RS           =   substructure od R&S parameters
%      .iz         =     indifference zone parameter
%      .alpha      =     alpha significance parameter
%*******************************************************************************
function results_gui(h,event)

% Retrieve application data and initialize
gui_func = getappdata(gcbf,'gui_func');
gui_var  = getappdata(gcbf,'gui_var');
gui      = getappdata(gcbf,'gui');
RunSet   = getappdata(gcbf,'RunSet');
k        = get(gcbo,'Position');
digits   = 10;
noYes    = {'No','Yes'};
BestF    = RunSet(k).BestF;
BestI    = RunSet(k).BestI;
isRS     = ~isempty(RunSet(k).RunStats.RS);

% Set up "Display Results" figure window
fig = figure( ...
   'Name',                           ['MADS Results for Problem: ', ...
                                     get(gui.problem, 'String'), ...
                                     ', Run #' int2str(k)], ...
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

% Set up panels on the figure window
if gui_var.newMatlab
   uipanel('Parent',fig,'Position',[0,  0,  1,.63], ...
                        'BorderType','beveledin','BackgroundColor',[.8 .8 .8]);
   uipanel('Parent',fig,'Position',[0, .63,.5,.37], ...
                        'BorderType','beveledin','BackgroundColor',[.8 .8 .8]);
   uipanel('Parent',fig,'Position',[.5,.63,.5,.37], ...
                        'BorderType','beveledin','BackgroundColor',[.8 .8 .8]);
else
   uicontrol(fig, 'Style','frame','Position', [0, .626, 1, .004],...
      'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
   uicontrol(fig, 'Style','frame','Position', [0, .63, 1, .004],...
      'ForegroundColor','white','BackgroundColor','white');
   uicontrol(fig, 'Style','frame','Position', [.5, .63, .003, .37], ...
      'ForegroundColor',[.4 .4 .4],'BackgroundColor',[.4 .4 .4]);
   uicontrol(fig, 'Style','frame','Position', [.497, .63, .003, .37], ...
      'ForegroundColor','white','BackgroundColor','white');
end

% Display Best Feasible Solution with pushbutton to view the solution
uicontrol(fig, 'String','BEST FEASIBLE SOLUTION', ...
         'FontWeight','bold', ...
         'HorizontalAlignment','center','Position',[.02 .86 .46 .10]);

if (isempty(RunSet(k).BestF))
   uicontrol(fig, 'String','NO FEASIBLE ITERATES FOUND', ...
         'HorizontalAlignment','center','Position',[.02 .77 .46 .08]);
else
   uicontrol(fig, 'String','Objective Function Value: ', ...
         'HorizontalAlignment','left',  'Position',[.02 .80 .28 .08]);
   uicontrol(fig, 'String','Constraint Violation Measure: ', ...
         'HorizontalAlignment','left',  'Position',[.02 .74 .28 .08]);
   uicontrol(fig, 'String',num2str(RunSet(k).BestF.f,digits), ...
         'HorizontalAlignment','right', 'Position',[.30 .80 .18 .08]);
   uicontrol(fig, 'String',num2str(RunSet(k).BestF.h,digits), ...
         'HorizontalAlignment','right', 'Position',[.30 .74 .18 .08]);
   uicontrol(fig,   'String',             'View Solution', ...
                    'Style',              'pushbutton', ...
                    'FontWeight',         'bold', ...
                    'HorizontalAlignment','center', ...
                    'Position',           [.02 .68 .46 .08], ...
                    'UserData',           'Best Feasible', ...
                    'Callback',           {gui_func.solution_gui,k,BestF});
end

% Display Least Infeasible Solution with pushbutton to view the solution
uicontrol(fig, 'String','LEAST INFEASIBLE SOLUTION', ...
         'FontWeight','bold', ...
         'HorizontalAlignment','center','Position',[.52 .86 .46 .10]);

if (isempty(RunSet(k).BestI))
   uicontrol(fig, 'String','NO INFEASIBLE ITERATES FOUND ', ...
         'HorizontalAlignment','center','Position',[.52 .77 .46 .08]);
else
   uicontrol(fig,'String','Objective Function Value: ', ...
         'HorizontalAlignment','left',  'Position',[.52 .80 .28 .08]);
   uicontrol(fig,'String','Constraint Violation Measure: ', ...
         'HorizontalAlignment','left',  'Position',[.52 .74 .28 .08]);
   uicontrol(fig,'String',num2str(RunSet(k).BestI.f,digits), ...
         'HorizontalAlignment','right', 'Position',[.80 .80 .18 .08]);
   uicontrol(fig,'String',num2str(RunSet(k).BestI.h,digits), ...
         'HorizontalAlignment','right', 'Position',[.80 .74 .18 .08]);
   uicontrol(fig,   'String',             'View Solution', ...
                    'Style',              'pushbutton', ...
                    'FontWeight',         'bold', ...
                    'HorizontalAlignment','center', ...
                    'Position',           [.52 .68 .46 .08], ...
                    'UserData',           'Least Infeasible', ...
                    'Callback',           {gui_func.solution_gui,k,BestI});
end

% Display Labels for MADS Run Statistics
uicontrol(fig, 'String','RUN STATISTICS', ...
   'FontWeight','bold', ...
   'HorizontalAlignment','center','Position',[.30 .51 .40 .10]);
uicontrol(fig, 'String','Final Mesh Size:', ...
   'HorizontalAlignment','left',  'Position',[.20 .45 .35 .08]);
uicontrol(fig, 'String','MADS Iterations:', ...
   'HorizontalAlignment','left',  'Position',[.20 .40 .35 .08]);
uicontrol(fig, 'String','Function Evaluations:', ...
   'HorizontalAlignment','left',  'Position',[.20 .35 .35 .08]);
uicontrol(fig, 'String','Gradient Evaluations:', ...
   'HorizontalAlignment','left',  'Position',[.20 .30 .35 .08]);
uicontrol(fig, 'String','CPU Time:', ...
   'HorizontalAlignment','left',  'Position',[.20 .25 .35 .08]);
uicontrol(fig, 'String','Consecutive Poll Failures:', ...
   'HorizontalAlignment','left',  'Position',[.20 .20 .35 .08]);
uicontrol(fig, 'String','Interrupted by User:', ...
   'HorizontalAlignment','left',  'Position',[.20 .15 .35 .08]);
uicontrol(fig, 'String','Cache Hits:', ...
   'HorizontalAlignment','left',  'Position',[.20 .10 .35 .08]);
if isRS
   uicontrol(fig, 'String','R&S Indifference Zone Parameter:', ...
      'HorizontalAlignment','left',  'Position',[.20 .05 .35 .08]);
   uicontrol(fig, 'String','R&S Alpha Parameter:', ...
      'HorizontalAlignment','left',  'Position',[.20 .00 .35 .08]);
end

% Display MADS Run Statistics (Compare with termination criteria)
Run.delta       = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.delta,digits), ...
   'HorizontalAlignment','right',  'Position',[.55 .45 .20 .08]);
Run.nIter       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.nIter), ...
   'HorizontalAlignment','right',  'Position',[.55 .40 .20 .08]);
Run.nFunc       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.nFunc), ...
   'HorizontalAlignment','right',  'Position',[.55 .35 .20 .08]);
Run.nGrad       = uicontrol(fig, ...
   'String',int2str(RunSet(k).RunStats.grad), ...
   'HorizontalAlignment','right',  'Position',[.55 .30 .20 .08]);
Run.time        = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.time), ...
   'HorizontalAlignment','right',  'Position',[.55 .25 .20 .08]);
Run.nFails      = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.nFails), ...
   'HorizontalAlignment','right',  'Position',[.55 .20 .20 .08]);
Run.stop        = uicontrol(fig, ...
   'String',noYes{1+RunSet(k).RunStats.stopRun}, ...
   'HorizontalAlignment','right',  'Position',[.55 .15 .20 .08]);
Run.nCacheHits  = uicontrol(fig, ...
   'String',num2str(RunSet(k).RunStats.nCacheHits), ...
   'HorizontalAlignment','right',  'Position',[.55 .10 .20 .08]);
if isRS
   Run.RS.iz     = uicontrol(fig, ...
      'String',num2str(RunSet(k).RunStats.RS.iz), ...
      'HorizontalAlignment','right',  'Position',[.55 .05 .20 .08]);
   Run.RS.alpha  = uicontrol(fig, ...
      'String',num2str(RunSet(k).RunStats.RS.alpha), ...
      'HorizontalAlignment','right',  'Position',[.55 .00 .20 .08]);
end

% Change colors for violated Termination Criteria
if (RunSet(k).RunStats.delta < str2double(get(gui.Term.delta,'String')))
   set(Run.delta,'FontWeight','bold','ForegroundColor','blue');
end
if (RunSet(k).RunStats.nIter >= str2double(get(gui.Term.nIter,'String')))
   set(Run.nIter, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.nFunc >= str2double(get(gui.Term.nFunc,'String')))
   set(Run.nFunc, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.time >= str2double(get(gui.Term.time,'String')))
   set(Run.time, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.nFails >= str2double(get(gui.Term.nFails,'String')))
   set(Run.nFails, 'FontWeight','bold','ForegroundColor','red');
end
if (RunSet(k).RunStats.stopRun)
   set(Run.stop, 'FontWeight','bold','ForegroundColor','red');
end
return

%*******************************************************************************
% solution_gui: View MADS optimal solution.
% ------------------------------------------------------------------------------
% Called by: results_gui (callback: View Solution pushbutton)
% Calls:     viewPrevious, viewNext
% VARIABLES:
%  iterate     = the best solution found by MADS
%    .x        =   vector of continuous variable values
%    .p        =   cell array of categorical variables
%  fig         = handle for the NOMADm GUI
%  gui_func    = handles for the NOMADm functions
%  gui_var     = handles for all GUI variables
%  gui.problem = name of current optimization problem
%  label       = text for "Best Feasible" or "Least Infeasible" solution
%  nX          = number of continuous variables
%  nP          = number of categorical variables
%  maxPage     = maximum number of variables displayed per page
%  nPages      = number of pages to display
%  onoff       = cell array containing "on" and "off"
%  leftmargin  = left margin position of continuous variables display
%  catmargin   = left margin position of categorical variables display
%  pos         = vertical position of next variable
%  xDisplay    = string conversion of the continuous variables
%  x(k)        = handles of the k-th displayed continuous variable
%    .ind      =   handle for the index label
%    .val      =   handle for the variable values
%  previousX   = handle for the continuous variables Previous push button
%  nextX       = handle for the continuous variables Next push button
%  pDisplay    = string conversion of the categorical variables
%  p(k)        = handles of the k-th displayed categorical variable
%    .ind      =   handle for the index label
%    .val      =   handle for the variable values
%  previousP   = handle for the categorical variables Previous push button
%  nextP       = handle for the categorical variables Next push button
%*******************************************************************************
function solution_gui(h,event,k,iterate)

fig      = findobj('Tag', 'NOMADm_GUI');
gui_func = getappdata(fig,'gui_func');
gui_var  = getappdata(fig,'gui_var');
gui      = getappdata(fig,'gui');
label    = get(gcbo,'UserData');
nX       = length(iterate.x);
nP       = length(iterate.p);
maxPage  = 12;
onoff    = {'on','off'};

% Set up "Display Solution" figure window
fig = figure('Name', ['MADS ',label,' Solution for Problem: ', ...
                      get(gui.problem, 'String'),', Run #' int2str(k)], ...
   'DefaultUIControlUnits',          'normalized', ...
   'DefaultUIControlFontUnits',      'normalized', ...
   'DefaultUIControlFontName',       'Helvetica',  ...
   'DefaultUIControlFontSize',        0.375, ...,
   'DefaultUIControlStyle',           'text', ...
   'Units',                           'normalized', ...
   'Position',                        [0.15 0.10 0.8 0.7], ...
   'MenuBar',                         'none', ...
   'NumberTitle',                     'off');

if gui_var.newMatlab
   uipanel;
else
   panelColor = get(fig,'DefaultUIControlBackgroundColor');
   uicontrol(fig, 'Style','frame','Position',[0, 0, 1, 1], ...
      'ForegroundColor',panelColor,'BackgroundColor',panelColor);
end

% Set display margins
if nP
   leftmargin = .10;
   catmargin  = .58;
else
   leftmargin = .30;
   catmargin  = .58;
end

% Display continuous variables Solution Labels
uicontrol(fig, 'String',             'CONTINUOUS VARIABLES', ...
               'FontWeight',         'bold', ...
               'HorizontalAlignment','center', ...
               'Position',           [leftmargin .86 .32 .10]);
uicontrol(fig, 'String',             'Index', ...
               'HorizontalAlignment','left', ...
               'Position',           [leftmargin .76 .07 .10]);
uicontrol(fig, 'String',             'Value', ...
               'HorizontalAlignment','right', ...
               'Position',           [leftmargin+.09 .76 .23 .10]);

% Display continuous variables
[xDisplay{1:maxPage}] = deal('');
for k = 1:nX
   xDisplay{k} = num2str(iterate.x(k));
end

for k = 1:maxPage
   pos = .76 - .05*k;
   x(k).ind = uicontrol(fig, ...
                        'String',             int2str(k), ...
                        'HorizontalAlignment','left', ...
                        'Position',           [leftmargin pos .05 .08]);
   x(k).val = uicontrol(fig, ...
                        'String',             xDisplay{k}, ...
                        'HorizontalAlignment','right', ...
                        'Position',           [leftmargin+.07 pos .25 .08]);
end
set([x(nX+1:end).ind, x(nX+1:end).val]', 'Visible','off');

% Display Previous/Next pushbuttons for continuous variables
previousX = uicontrol(fig, ...
                      'Style',    'pushbutton', ...
                      'String',   'Previous',   ...
                      'Visible',  'off', ...
                      'Position', [leftmargin+.07, .05, .08, .06], ...
                      'Callback', {gui_func.changeView});
nextX     = uicontrol(fig, ...
                      'Style',    'pushbutton', ...
                      'String',   'Next', ...
                      'Visible',  onoff{2 - (nX > maxPage)}, ...
                      'Position', [leftmargin+.19, .05, .08, .06], ...
                      'Callback', {gui_func.changeView});
set(previousX,'UserData',{-1,x,iterate.x,nextX});
set(nextX,    'UserData',{ 1,x,iterate.x,previousX});

% If MVP problem:
if nP

   % Convert the categorical variables into a cell array of strings
   [pDisplay{1:maxPage}] = deal('');
   for k = 1:nP
      if ischar(iterate.p{k})
         pDisplay{k} = iterate.p{k};
      else
         pDisplay{k} = num2str(iterate.p{k});
      end
   end

   % Display categorical variables Solution Labels
   uicontrol(fig, 'String',             'CATEGORICAL VARIABLES', ...
                  'FontWeight',         'bold', ...
                  'HorizontalAlignment','center', ...
                  'Position',           [catmargin .86 .32 .10]);
   uicontrol(fig, 'String',             'Index', ...
                  'HorizontalAlignment','left', ...
                  'Position',           [catmargin .76 .07 .10]);
   uicontrol(fig, 'String',             'Value', ...
                  'HorizontalAlignment','right', ...
                  'Position',           [catmargin+.09 .76 .23 .10]);

   % Display categorical variables
   for k = 1:maxPage
      pos = .76 - .05*k;
      p(k).ind = uicontrol(fig, ...
                           'String',             int2str(k), ...
                           'HorizontalAlignment','left', ...
                           'Position',           [catmargin pos .05 .08]);
      p(k).val = uicontrol(fig, ...
                           'String',             pDisplay{k}, ...
                           'HorizontalAlignment','right', ...
                           'Position',           [catmargin+.07 pos .25 .08]);
   end
   set([p(nP+1:end).ind, p(nP+1:end).val]', 'Visible','off');

   % Display Previous/Next pushbuttons for categorical variables
   previousP = uicontrol(fig, ...
                         'Style',    'pushbutton', ...
                         'String',   'Previous',   ...
                         'Visible',  'off', ...
                         'Position', [catmargin+.07, .05, .08, .06], ...
                         'Callback', {gui_func.changeView});
   nextP     = uicontrol(fig, ...
                         'Style',    'pushbutton', ...
                         'String',   'Next', ...
                         'Visible',  onoff{2 - (nP > maxPage)}, ...
                         'Position', [catmargin+.19, .05, .08, .06], ...
                         'Callback', {gui_func.changeView});
   set(previousP,'UserData',{-1,p,iterate.p,nextP});
   set(nextP,    'UserData',{ 1,p,iterate.p,previousP});
end
return

%*******************************************************************************
% changeView: Change the view of a multi-page list of variable values.
% ------------------------------------------------------------------------------
% Called by: solutions_gui (callback: Previous/Next pushbuttons)
% VARIABLES:
%   maxPage     = maximum number of variables displayed on one page
%   UserData    = storage for key items needed
%   code        = an integer: magnitude = page #, sign = which button pressed
%   solution    = handles for the on-screen display of the variable values
%     .ind      =   handles for the variable indices
%     .val      =   handles for the variable values
%   iterate     = the solution vector
%   otherhandle = handle of the push button not pressed
%   whichButton = pushbutton indicator: Next = +1, Previous = -1
%   iPage       = display page number
%   n           = number of variables
%   lastPage    = the total number of display pages
%   newcode     = updated code variable assigned back to the pushbutton UserData
%   nDisplay    = number of variables displayed on the current page
%   ind         = index of variable to be displayed
%   value       = variable value converted to string for display
%*******************************************************************************
function changeView(h,event)

% Retrieve pushbutton parameters
UserData = get(gcbo, 'UserData');
[code,solution,iterate,otherhandle] = deal(UserData{:});

% Set up page parameters and update pushbutton parameters for new page
whichButton = sign(code);
iPage       = abs(code);
iPage       = iPage + whichButton;
maxPage     = length(solution);
n           = length(iterate);
lastPage    = ceil(n/maxPage);
newcode     = iPage*whichButton;
set(gcbo,        'UserData',{ newcode,solution,iterate,otherhandle});
set(otherhandle, 'UserData',{-newcode,solution,iterate,gcbo});

% Set visibility of variable and index displays
set([gcbo; otherhandle], 'Visible', 'on');
nDisplay = maxPage;
if (iPage == 1)
   set(gcbo, 'Visible','off');
end
if (iPage == lastPage)
   nDisplay = max(1,mod(n,maxPage));
   set(gcbo, 'Visible','off');
end

% Update display of indices and variables for new page
for k = 1:nDisplay
   ind = k + (iPage-1)*maxPage;
   if iscell(iterate)
      if ischar(iterate{ind})
         value = iterate{ind};
      else
         value = num2str(iterate{ind});
      end
   else
      value = num2str(iterate(ind));
   end
   set(solution(k).ind, 'Visible','on','String',int2str(ind));
   set(solution(k).val, 'Visible','on','String',value);
end
for k = nDisplay+1:maxPage
   set([solution(k).ind; solution(k).val], 'Visible','off');
end
return
