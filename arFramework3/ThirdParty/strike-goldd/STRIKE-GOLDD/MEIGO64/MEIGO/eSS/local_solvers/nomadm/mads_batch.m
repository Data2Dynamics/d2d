function [BestF,BestI,RunStats,RunSet] = mads_batch(problem_path,tol_nomad,nvar,nconst,iterprint)
%MADS_BATCH  Sets up and runs the MADS algorithm without a GUI.
%
%   Syntax:
%      mads_batch
%
%   Description:
%      This function serves as a GUI-free alternative to NOMADm in setting
%      up an optimization problem, setting various algorithm parameters and
%      user options, and calling the MADS optimizer.  It first sets all of the
%      variables to their default values, which are clearly stated in the
%      MADS_DEFAULTS file.  To change a variable from its default value, the
%      user must add a statement to this file to do so.  Some variable
%      statements are included here for convenience, which can be change
%      manually.
%
%   See also MADS_DEFAULTS, MADS

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
%   You should have received a copy of the GNU General Public License along
%   with NOMADm; if not, write to the Free Software Foundation, Inc., 
%   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% ------------------------------------------------------------------------------
%   Originally created, 2001.
%   Last modified, 31 January 2005
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
% mads_batch:  Runs MADS in batch mode.
% ------------------------------------------------------------------------------
% Calls:     mads_defaults, mads, < user initial points file >
% Variables:
%  Defaults    = structure of MADS default values (see mads_defaults)
%  Options     = structure for options settings (see mads_defaults)
%  problemPath = location of user problem files
%  Problem     = structure of data for optimization problem
%  newPath     = logical indicating if path is not the Matlab path
%  iterate0    = structure of data for the initial iterate (see mads)
%  BestF       = final best feasible solution found
%  BestI       = final least infeasible solution found
%  RunStats    = structure of MADS Run statistics (see mads)
%*******************************************************************************

% Set Options to their default values
%clear variables
Defaults = mads_defaults('Truth');
Options  = Defaults.Options;
Problem.nameCache = 'CACHE';

%*******************************************************************************
% MODIFY ONLY AFTER THIS POINT
%*******************************************************************************

% Specify Problem Files
problemPath=problem_path;


Problem.File.F = 'fobj_nomad';           % functions file
Problem.File.O = 'fobj_nomad_Omega';     % Omega file (linear constraints)
Problem.File.I = 'fobj_nomad_x0';        % initial Points file
Problem.File.N = 'fobj_nomad_N';         % discrete neighbor file (MVP only)
Problem.File.C = 'fobj_nomad_cache.mat'; % previously created Cache file
Problem.File.P = 'fobj_nomad_Param';     % previously created Session file
Problem.fType  = 'M';                   % type of Functions file {M,F,C}
Problem.nc     = 0;                     % number of nonlinear constraints


% Specify Choices for SEARCH
Options.Search(1).type    = 'None';      % For choices, see mads_defaults
Options.Search(1).nIter   = 0;           % Number of iterations for Search #1
Options.Search(1).nPoints = 0;           % Number of poll or sample points
Options.Search(1).file    = '';          % filename must include the full path
Options.Search(1).local   = 0;           % flag to turn on trust region
Options.Search(1).merit   = 0;           % flag to penalize clustered data
Options.Search(2).type    = 'None';      % For choices, see mads_defaults
Options.Search(2).nIter   = Inf;         % Number of iterations for Search #2
Options.Search(2).nPoints = 0;           % Number of poll or sample points
Options.Search(2).file    = '';          % filename must include the full path
Options.Search(2).local   = 0;           % flag to turn on trust region
Options.Search(2).merit   = 0;           % flag to penalize clustered data
Options.nSearches = length(Options.Search);
Options.SurOptimizer      = 'fmincon';

% Specify Choices for POLL
Options.pollStrategy    = 'Standard_2n'; % For choices, see mads_defaults
Options.pollOrder       = 'Consecutive'; % For choices, see mads_defaults
Options.pollCenter      = 0;             % Poll around n-th filter point
Options.pollComplete    = 0;             % Flag for complete polling

% Specify Termination Criteria

%==========================================================================
%JAE 09/05/07
Options.Term.delta      = tol_nomad;      %mesh size
%==========================================================================
Options.Term.nIter      = Inf;           % maximum number of iterations
Options.Term.nFunc      = 100*nvar;         % maximum number of function evals
Options.Term.time       = Inf;           % maximum CPU time
Options.Term.nFails     = Inf;           % max number of consecutive Poll fails

% Choices for Mesh Control
Options.delta0          = 1.0;           % initial mesh size
Options.deltaMax        = Inf;           % bound on how coarse the mesh can get
Options.meshRefine      = 0.5;           % mesh refinement factor
Options.meshCoarsen     = 2.0;           % mesh coarsening factor

% Choices for Filter management (for problems with nonlinear constraints)
if ~nconst, nconst=1; end
Options.hmin            = nconst*tol_nomad*tol_nomad;          % minimum infeasible point h-value
%Options.hmin            =  1e-5;          % minimum
%infeasible point h-value
Options.hmax            = 1.0;           % maximum h-value of a filter point

% Choices for EXTENDED POLL (for MVP problems)
Options.ePollTriggerF   = 0.01;          % f-value Extended Poll trigger
Options.ePollTriggerH   = 0.01;          % h-value Extended Poll trigger

% MADS flag parameter values
Options.loadCache        = 1;            % load pre-existing Cache file
Options.countCache       = 1;            % count Cache points as function calls
Options.runStochastic    = 0;            % runs problem as a stochastic problem
Options.removeRedundancy = 1;            % discard redundant linear constraints
Options.scale            = 2;            % scale directions using this log base
Options.computeGrad      = 0;            % compute gradient, if available
Options.plotFilter       = 0;            % plot the filter real-time

%==========================================================================
%JAE 26/08/05
if iterprint
    Options.plotHistory1     = 1;            % plot MADS performance 
    Options.plotHistory2     = 1;            % plot MADS performance real-time
else
    Options.plotHistory1     = 0;            % plot MADS performance 
    Options.plotHistory2     = 0;            % plot MADS performance real-time
end
%==========================================================================
    
Options.plotColor        = 'k';          % color of history plot

%*******************************************************************************
% DO NOT MODIFY AFTER THIS POINT
%*******************************************************************************

% Set up figure handles for real-time plots
if (Problem.nc == 0)
   Options.plotFilter = 0;
end
if (Options.plotFilter)
   figure;
   Options.fplothandle = gca;
end
if (Options.plotHistory2)
   figure;
   Options.hplothandle = gca;
end
% Set the path, get the initial iterates, call the optimizer
newPath = isempty(strfind(upper(path),upper(problemPath)));
addpath(problemPath);
iterate0 = feval(Problem.File.I);
[BestF,BestI,RunStats,RunSet] = mads(Problem,iterate0,Options);
if newPath, rmpath(problemPath); end
return
