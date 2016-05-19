% This function compiles the Ceres code with the mex wrapper.
% It should be called to produce a mex file which is callable. Note that this requires
% a C compiler


%% LICENSE MESSAGE FROM CERES %%
% Ceres Solver - A fast non-linear least squares minimizer
% Copyright 2015 Google Inc. All rights reserved.
% http://ceres-solver.org/
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice,
%   this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% * Neither the name of Google Inc. nor the names of its contributors may be
%   used to endorse or promote products derived from this software without
%   specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%% END LICENSE MESSAGE FROM CERES %%


function compileCeres ()   
    % Check if we have a C compiler
    try
        cc = mex.getCompilerConfigurations('C++');
        fprintf( 'C++ compiler(s):\t\t%s\n', sprintf('[%s] ', cc.Name) );
    catch
        error( 'No C++ Compiler installed' );
    end    
    
    cpath   = mfilename('fullpath');
    if (ispc)
        slash = '\';
        objectExtension = '.obj';
    else
        slash = '/';
        objectExtension = '.o';
    end
    loc     = strfind( fliplr(cpath), slash);
    cpath   = cpath(1:end-loc+1);
    
%% %%%%%%%%%%%%%%%%%%%%  CERES FILE LIST OF CC FILES %%%%%%%%%%%%%%%%%%%%%
    
    fileListCeres    = getAllFiles(strcat(sprintf('%s',cpath),[slash 'ceres-solver']));
    
    preselect        = strfind(fileListCeres,'.cc');
    pre1select       = strfind(fileListCeres, 'test_util.cc');
    pre2select       = strfind(fileListCeres, 'collections_port.cc');
    pre3select       = strfind(fileListCeres, 'gmock_main.cc');
    pre4select       = strfind(fileListCeres, 'autodiff');
    pre5select       = strfind(fileListCeres, 'numeric_diff');    
    pre6select       = strfind(fileListCeres, 'matrix_utils_test');
    pre7select       = strfind(fileListCeres, 'compressed_row_sparse_matrix_test');    
    pre8select       = strfind(fileListCeres, 'dense_sparse_matrix_test');
    pre9select       = strfind(fileListCeres, 'jet_test');
    pre10select      = strfind(fileListCeres, 'bundle_adjustment_test');    
    pre11select      = strfind(fileListCeres, 'polynomial_test');
    pre12select      = strfind(fileListCeres, 'rotation_test');
    
    
    select = (~isnan(cellfun(@mean,preselect)));
    
    p1 = isnan(cellfun(@mean,pre1select));
    p2 = isnan(cellfun(@mean,pre2select));
    p3 = isnan(cellfun(@mean,pre3select));
    p4 = isnan(cellfun(@mean,pre4select));
    p5 = isnan(cellfun(@mean,pre5select));
    p6 = isnan(cellfun(@mean,pre6select));
    p7 = isnan(cellfun(@mean,pre7select));
    p8 = isnan(cellfun(@mean,pre8select)); 
    p9 = isnan(cellfun(@mean,pre9select));
    p10 = isnan(cellfun(@mean,pre10select));
    p11 = isnan(cellfun(@mean,pre11select));
    p12 = isnan(cellfun(@mean,pre12select));
    
    ff = logical(p1.*p2.*p3.*p4.*p5.*p6.*p7.*p8.*p9.*p10.*p11.*p12.*select); 
    
    ccfileListCeres = fileListCeres(ff);
    
    ccoutFilesCeres = ccfileListCeres;
    
    
    
    for i = 1:length(ccoutFilesCeres)
            loc             = strfind( fliplr(ccoutFilesCeres{i}), slash);
            ccoutFilesCeres{i}   = ccoutFilesCeres{i}(end-loc+2:end);  
            ccoutFilesCeres{i}   = strrep(ccoutFilesCeres{i},'.cc',objectExtension);
    end
   
    %ccoutFiles = cellfun(@strrep('.cc', '.o'), ccoutFiles);
    
   %%%%%%%%%%%%%%%%%%%%  END CERES FILE LIST OF CC FILES %%%%%%%%%%%%%%%%%%
   
   
    
%% %%%%%%%%%%%%%%%%%%%% SuiteSparse CHOLMOD %%%%%%%%%%%%%%%%%%%
    
%     fileListCHOLMOD = getAllFiles(strcat(sprintf('%s',cpath),'/SuiteSparse/CHOLMOD/Core'));
%     preselectcholmod       = strfind(fileListCHOLMOD,'.c');
%     selectCHOLMOD = (~isnan(cellfun(@mean,preselectcholmod)));
%     ccfileListCHOLMOD = fileListCHOLMOD(selectCHOLMOD);
%     
%     ccoutFilesCHOLMOD = ccfileListCHOLMOD;
% 
%     
%     for i = 1:length(ccoutFilesCHOLMOD)
%             loc             = strfind( fliplr(ccoutFilesCHOLMOD{i}), '/');
%             ccoutFilesCHOLMOD{i}   = ccoutFilesCHOLMOD{i}(end-loc+2:end);  
%             ccoutFilesCHOLMOD{i}   = strrep(ccoutFilesCHOLMOD{i},'.cc','.o');
%     end

% %%%%%%%%%%%%%%%%%%%% SuiteSparse CHOLMOD END %%%%%%%%%%%%%%%%%%%    
    
    


%% Compilation

    fprintf( 'Compiling Ceres... \t' );
    
    includesstr = {};
    
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinternal%sceres"', slash, slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinternal%sceres%sminiglog"', slash, slash, slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinternal"', slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinclude"', slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%seigen3"', slash)];
    
    disp('compiling');
    mex  ( '-c', 'CXXFLAGS=''\$CXXFLAGS -fPIC''', includesstr{:}, '-largeArrayDims', '-lmwblas', '-lmwlapack', ccfileListCeres{:}  );
    disp('linking');
    mex  ( includesstr{:}, 'CXXFLAGS=''\$CXXFLAGS -fPIC''', '-largeArrayDims', '-lmwblas', '-lmwlapack',sprintf('%s%sceresd2d.cpp', cpath, slash), ccoutFilesCeres{:}, '-outdir', cpath);
    
%    mex  -I"/usr/include/eigen3"  -largeArrayDims -lmwblas -lstdc++ -L"/home/fgwieland/Bachelor/Matlab-Arbeiten/Eigene Matlab Arbeiten/Solver/CERES/Eigene Implementierung/d2dImplementierung" -lceresd2d '/home/fgwieland/Bachelor/Matlab-Arbeiten/Eigene Matlab Arbeiten/Solver/CERES/Eigene Implementierung/d2dImplementierung/ceresd2d.cpp'
%    mex  -I"/usr/include/eigen3"  -largeArrayDims -lmwblas -L"/usr/lib/x86_64-linux-gnu/" -lstdc++ -L"/home/fgwieland/Bachelor/Matlab-Arbeiten/Eigene Matlab Arbeiten/Solver/CERES/Eigene Implementierung/d2dImplementierung" -lceresd2d '/home/fgwieland/Bachelor/Matlab-Arbeiten/Eigene Matlab Arbeiten/Solver/CERES/Eigene Implementierung/d2dImplementierung/ceresd2d.cpp'


%%  Cleanup and exit message
   %  delete(ccoutFilesCeres{:});

    fprintf( '[ CERES successfully compiled ]\n');    
end


%% Function to get all files from folder & subfolders
function fileList = getAllFiles(dirName)

  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
  end

end


