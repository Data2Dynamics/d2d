% NOMADm Software Package
%   Version 3.3, February 2, 2005
%   Copyright (c) 2001-2005 by Mark A. Abramson
%
%   NOMADm files:
%      nomadm          - Executes the graphical user interface (GUI) and
%                        callback functions for running the MADS optimizer.
%      mads_batch      - Sets up user options and runs MADS without the GUI.
%      mads            - Runs the mesh-adaptive direct search (MADS) filter
%                        algorithm for solving constrained or unconstrained
%                        nonlinear or mixed variable optimization problems.
%      mads_defaults   - Sets the default values used by all the NOMADm files.
%      nomadm_compile  - A utility that uses mex to compile a C/C++ or FORTRAN
%                        user functions file so it can be called in MATLAB.
%
%   Documentation:
%      nomadm_help.pdf - NOMADm User's Guide.
%      nomadm.txt      - Summary of all software changes in recent versions.
%      gpl.txt         - GNU General Public License.
%
%   Other packages (URLs provided below):
%      dace            - DACE Toolbox - used for forming surrogate problems
%      nw              - a package for constructing surrogates using the
%                        Nadaraya-Watson estimator.
%      cmaes           - Genetic algorithm - used to optimize surrogates
%                        and as a Search option
%
%      dace:  http://www.imm.dtu.dk/~hbn/dace/
%      nw:    http://en.afit.edu/enc/Faculty/MAbramson/nomadm.html
%      cmaes: http://www.bionik.tu-berlin.de/user/niko/cmaes_inmatlab.html

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
%   Last modified, 28 January 2005
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
