function language = nomadm_compile(fileType,filePath,fileName,fileExt)
%NOMADM_COMPILE  Uses MEX to compile a C and FORTRAN functions file.
%
%   Syntax:
%      LANGUAGE = nomadm_compile(FILETYPE,FILEPATH,FILENAME,FILEEXT)
%
%   Desciption:
%      NOMADM_COMPILE is a utility for compiling functions files that are
%      written in C/C++ or FORTRAN for use with the MADS optimizer.  This can
%      be done manually by following the instructions here, or automatically
%      through the NOMADm GUI, in which case, it is called by NOMADM.  LANGUAGE
%      is a string used by the NOMADm software to identify the programming
%      language of the source code. If run separately, LANGUAGE is ignored.
%
%      FILETYPE is a 1-character string code that identifies the programming
%      language in which the functions file is written.  The legal codes are
%      F = FORTRAN and C = C/C++.  If called from nomadm, this code is obtained
%      by taking the first character of the filename extension.  The variables,
%      FILEPATH, FILENAME, and FILEEXT help to identify the location, filename,
%      and filename extension, respectively, of the functions file.
%
%   See also MADS, NOMADM

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
% nomadm_compile: Construct a MEX gateway file for a non-MATLAB function file.
% ------------------------------------------------------------------------------
% Called by:  loadMADS
% Calls:      compileFortran, compileC
% VARIABLES:
%  language = programming language of Functions file
%  fileType = string code for functions file programming language
%  filePath = current path where the non-Matlab functions file lies
%  fileName = optimization problem name
%  fileExt  = filename extention of the non-Matlab functions file
%*******************************************************************************
switch upper(fileType)
   case 'F'
      compileFortran(filePath,fileName,fileExt);
      language = 'FORTRAN';
   case 'C'
      compileC(filePath,fileName,fileExt);
      language = 'C';
   case 'M'
      warning('Matlab Function file found.  No compiling is necessary.');
      language = 'Matlab';
   otherwise
      error('Unable to compile:  Incompatible File Type');
end
return

%*******************************************************************************
% compileFortran: Construct a MEX Fortran gateway file.
% ------------------------------------------------------------------------------
% Called by:  nomadm_compile
% Calls:      getMexExt
% VARIABLES:
%  filePath    = current path where the Fortran functions file lies
%  fileName    = optimization problem name
%  fortranExt  = filename extention of the Fortran functions file
%  userFile    = entire filename of the Fortran functions file
%  wrapperFile = entire filename of gateway file
%  execFile    = entire filename of the compiled executable functions file
%  p           = cell array containing the Fortran commands to be executed
%  fid         = file ID -- handle for the file that is being written
%*******************************************************************************
function compileFortran(filePath,fileName,fortranExt)

% Construct the names of the gateway file and compiled file
userFile    = [filePath,fileName,fortranExt];
wrapperFile = [filePath,fileName,'g',fortranExt];
execFile    = [filePath,fileName,getMexExt];

% Construct a gateway file if one does not already exist
if (~exist(wrapperFile,'file'))
   p{1}     = 'C';
   p{end+1} = 'C Usage: [nc,fx,cx,gradfx,gradcx] = ExamplePF(x,p1,p2)';
   p{end+1} = 'C';
   p{end+1} = '       SUBROUTINE MEXFUNCTION(NLHS,PLHS,NRHS,PRHS)';
   p{end+1} = '       INTEGER NLHS, NRHS, PLHS(*), PRHS(*)';
   p{end+1} = 'C';
   p{end+1} = '       INTEGER MXCREATEFULL,MXGETPR,MXGETM,MXGETN';
   p{end+1} = 'C';
   p{end+1} = 'C KEEP THE ABOVE PORTION FOR USE IN ALL FORTRAN MEX FILES';
   p{end+1} = 'C -------------------------------------------------------';
   p{end+1} = 'C';
   p{end+1} = '       INTEGER NCMAX,NXMAX,NPMAX,NC,NX,NP1,NP2';
   p{end+1} = '       INTEGER PNC,PX,PP,PF,PC,PDF,PDC';
   p{end+1} = '       PARAMETER (NCMAX = 10,NPMAX = 10, NXMAX = 100)';
   p{end+1} = '       REAL*8  P1(NPMAX), X(NXMAX)';
   p{end+1} = '       REAL*8  F,C(NCMAX),DF(NXMAX),DC(NXMAX,NCMAX)';
   p{end+1} = '       CHARACTER P2(NPMAX)';
   p{end+1} = 'C';
   p{end+1} = '       IF (NRHS .NE. 3 .OR. NLHS .NE. 5) THEN';
   p{end+1} = '          CALL MEXERRMSGTXT(''Invalid number of arguments'')';
   p{end+1} = '       ENDIF';
   p{end+1} = 'C';
   p{end+1} = 'C Get input data and set up input data pointers';
   p{end+1} = '       NX      = MAX(MXGETM(PRHS(1)),MXGETN(PRHS(1)))';
   p{end+1} = '       NP1     = MAX(MXGETM(PRHS(2)),MXGETN(PRHS(2)))';
   p{end+1} = '       NP2     = MAX(MXGETM(PRHS(3)),MXGETN(PRHS(3)))';
   p{end+1} = '       PX      = MXGETPR(PRHS(1))';
   p{end+1} = '       PP      = MXGETPR(PRHS(2))';
   p{end+1} = 'C';
   p{end+1} = '       IF (NX .GT. NXMAX .OR. MAX(NP1,NP2) .GT. NPMAX) THEN';
   p{end+1} = '          CALL MEXERRMSGTXT(''Array sizes too small'')';
   p{end+1} = '       ENDIF';
   p{end+1} = 'C';
   p{end+1} = 'C';
   p{end+1} = 'C Allocate memory and set up pointers for output variables';
   p{end+1} = '       PLHS(1) = MXCREATEFULL(1,1,0)';
   p{end+1} = '       PLHS(2) = MXCREATEFULL(1,1,0)';
   p{end+1} = '       PLHS(3) = MXCREATEFULL(1,NCMAX,0)';
   p{end+1} = '       PLHS(4) = MXCREATEFULL(NX,1,0)';
   p{end+1} = '       PLHS(5) = MXCREATEFULL(NX,NCMAX,0)';
   p{end+1} = '       PNC     = MXGETPR(PLHS(1))';
   p{end+1} = '       PF      = MXGETPR(PLHS(2))';
   p{end+1} = '       PC      = MXGETPR(PLHS(3))';
   p{end+1} = '       PDF     = MXGETPR(PLHS(4))';
   p{end+1} = '       PDC     = MXGETPR(PLHS(5))';
   p{end+1} = 'C';
   p{end+1} = 'C Load data into Fortran arrays and call the main subroutine';
   p{end+1} = '       CALL MXCOPYPTRTOREAL8(PX,X,NX)';
   p{end+1} = '       CALL MXCOPYPTRTOREAL8(PP,P1,NP1)';
   p{end+1} = '       STATUS = MXGETSTRING(PRHS(3),P2,NP2)';
   p{end+1} = '       IF (STATUS .EQ. 1 .AND. NP2 .NE. 0) THEN';
   p{end+1} = '          CALL MEXERRMSGTXT(''Character array too small'')';
   p{end+1} = '       ENDIF';
   p{end+1} =['       CALL ',fileName,'(NP1,NP2,NX,P1,P2,X,NC,F,C,DF,DC)'];
   p{end+1} = 'C';
   p{end+1} = 'C Load output into Matlab structures';
   p{end+1} = '       CALL MXCOPYREAL8TOPTR(DFLOAT(NC),PNC,1)';
   p{end+1} = '       CALL MXCOPYREAL8TOPTR(F,PF,1)';
   p{end+1} = '       CALL MXCOPYREAL8TOPTR(C,PC,NC)';
   p{end+1} = '       CALL MXCOPYREAL8TOPTR(DF,PDF,NX)';
   p{end+1} = '       CALL MXCOPYREAL8TOPTR(DC,PDC,NC*NX)';
   p{end+1} = '       RETURN';
   p{end+1} = '       END';

   fid = fopen(wrapperFile,'w');
   fprintf(fid,'%s\n',p{:});
   fclose(fid);
end

% Try to compile the files, if a compiled file does not already exist
if (~exist(execFile,'file'))
   mex(userFile,wrapperFile,'-outdir',filePath);
end
return

%*******************************************************************************
% processC:  Construct a MEX C gateway file.
% ------------------------------------------------------------------------------
% Called by:  nomadm_functions
% Calls:      getMexExt
% VARIABLES:
%  filePath    = current path where the C functions files lies
%  fileName    = optimization problem name
%  cExt        = filename extention of the C functions file
%  userFile    = entire filename of the C functions file
%  wrapperFile = entire filename of gateway file
%  execFile    = entire filename of the compiled executable functions file
%  p           = cell array containing the C commands to be executed
%  fid         = file ID -- handle for the file that is being written
%*******************************************************************************
function compileC(filePath,fileName,cExt)

% Construct the names of the gateway file and compiled file
userFile    = [filePath,fileName,cExt];
wrapperFile = [filePath,fileName,'temp',cExt];
execFile    = [filePath,fileName,getMexExt];

% Construct a gateway file if one does not already exist
if (~exist(execFile,'file'))
   p{1}     = '#include "mex.h"';
   p{2}     = ' ';
   p{3}     = '/* Usage: [nc,fx,cx,gradfx,gradcx] = ExamplePC(x,p1,p2) */';
   p{4}     = ' ';
   fid = fopen(userFile,'r');
   while ~feof(fid)
      p{end+1} = fgetl(fid);
   end
   fclose(fid);
   p{end+1} = ' ';
   p{end+1} = 'void mexFunction(int nlhs, mxArray *plhs[],';
   p{end+1} = '                 int nrhs, const mxArray *prhs[])';
   p{end+1} = '{';
   p{end+1} = 'unsigned nx, np1, np2;';
   p{end+1} = 'const unsigned ncmax = 100;';
   p{end+1} = 'double   *p1, *x, *pnc, *pf, *pc, *pdf, *pdc;';
   p{end+1} = 'char     *p2;';
   p{end+1} = ' ';
   p{end+1} = '/* ERROR CHECK the number of input and output arguments */';
   p{end+1} = 'if (nrhs != 3) mexErrMsgTxt("Must have 3 Input Arguments");';
   p{end+1} = 'if (nlhs != 5) mexErrMsgTxt("Must have 5 Output Arguments");';
   p{end+1} = ' ';
   p{end+1} = '/* GET INPUT DATA and set up input data pointers */';
   p{end+1} = 'nx  = (unsigned) mxGetNumberOfElements(prhs[0]);';
   p{end+1} = 'np1 = (unsigned) mxGetNumberOfElements(prhs[1]);';
   p{end+1} = 'np2 = (unsigned) mxGetNumberOfElements(prhs[2]);';
   p{end+1} = 'x   = mxGetPr(prhs[0]);';
   p{end+1} = 'p1  = mxGetPr(prhs[1]);';
   p{end+1} = 'p2  = mxGetChars(prhs[2]);';
   p{end+1} = ' ';
   p{end+1} = '/* SET UP POINTERS AND MEMORY FOR OUTPUT VARIABLES */';
   p{end+1} = 'plhs[0] = mxCreateScalarDouble(0);';
   p{end+1} = 'plhs[1] = mxCreateScalarDouble(0);';
   p{end+1} = 'plhs[2] = mxCreateDoubleMatrix(1,ncmax,mxREAL);';
   p{end+1} = 'plhs[3] = mxCreateDoubleMatrix(nx,1,mxREAL);';
   p{end+1} = 'plhs[4] = mxCreateDoubleMatrix(nx,ncmax,mxREAL);';
   p{end+1} = 'pnc     = mxGetPr(plhs[0]);';
   p{end+1} = 'pf      = mxGetPr(plhs[1]);';
   p{end+1} = 'pc      = mxGetPr(plhs[2]);';
   p{end+1} = 'pdf     = mxGetPr(plhs[3]);';
   p{end+1} = 'pdc     = mxGetPr(plhs[4]);';
   p{end+1} = ' ';
   p{end+1} = '/* GET FUNCTION OUTPUT DATA */';
   p{end+1} = [fileName,'(np1,np2,nx,p1,p2,x,pnc,pf,pc,pdf,pdc);'];
   p{end+1} = '}';

   % Write the wrapper file source code
   fid = fopen(wrapperFile,'w');
   fprintf(fid,'%s\n',p{:});
   fclose(fid);

   % Compile the wrapper file and delete its source code
   mex(wrapperFile,'-outdir',filePath,'-output', fileName);
   if (exist(wrapperFile,'file')), delete(wrapperFile); end
end
return

%*******************************************************************************
% GetMexExt:  Returns platform-dependent extension for an executable file.
% ------------------------------------------------------------------------------
% Called by: compileFortran, compileC
% VARIABLES:
%  mexext   = platform-dependent executable file extension
%  computer = Matlab built-in string function indicating platform type
%*******************************************************************************
function mexext = getMexExt

% Determine the platform-dependent extention of the compiled file
switch computer
case {'PCWIN','MAC2'}
   mexext = '.dll';     % MS Windows, Mac
case {'LNX86'}
   mexext = '.mexglx';  % Linux Intel
case {'SUN4','SOL2'}
   mexext = '.mexsol';  % SUN Sparc, Solaris 2
case {'SGI','SGI64'}
   mexext = '.mexsg';   % Silicon Graphics
case {'IBM_RS'}
   mexext = '.mexrs6';  % IBM RS6000
case {'ALPHA','AXP_VMSG','AXP_VMSIEEE','VAX_VMSG','VAX_VMSD'}
   mexext = '.mexaxp';  % DEC Alpha, VAX/VMS
case {'HP700'}
   mexext = '.mexhp7';  % HP 9000/700 (also .mexhpx)
otherwise
   error('Unknown computer type');
end
return
