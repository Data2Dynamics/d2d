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
    
    select = (~isnan(cellfun(@mean,preselect)));
    
    excludelist = {};
    excludelist{end+1} = 'test_util.cc';
    excludelist{end+1} = 'collections_port.cc';
    excludelist{end+1} = 'gmock_main.cc';
    excludelist{end+1} = 'autodiff';
    excludelist{end+1} = 'numeric_diff';
    excludelist{end+1} = 'matrix_utils_test';
    excludelist{end+1} = 'compressed_row_sparse_matrix_test';
    excludelist{end+1} = 'dense_sparse_matrix_test';
    excludelist{end+1} = 'jet_test';
    excludelist{end+1} = 'bundle_adjustment_test';
    excludelist{end+1} = 'polynomial_test';
    excludelist{end+1} = 'rotation_test';
    
    ex = ones(length(fileListCeres),1);
    
    for i = 1:12
        temp    = strfind(fileListCeres, excludelist{i});
        ptemp   = isnan(cellfun(@mean,temp));
        
        ex      = ex.*ptemp;
    end
    
    ff = logical(ex.*select);

    ccfileListCeres = fileListCeres(ff);
    
    ccoutFilesCeres = ccfileListCeres;
    
    
    
    for i = 1:length(ccoutFilesCeres)
            loc             = strfind( fliplr(ccoutFilesCeres{i}), slash);
            ccoutFilesCeres{i}   = ccoutFilesCeres{i}(end-loc+2:end);  
            ccoutFilesCeres{i}   = strrep(ccoutFilesCeres{i},'.cc',objectExtension);
    end
   %%%%%%%%%%%%%%%%%%%%  END CERES FILE LIST OF CC FILES %%%%%%%%%%%%%%%%%%
 
%% Compilation
    fprintf( 'Copyright 2015 Google Inc. All rights reserved.\n');
    fprintf( 'Compiling Ceres... \t' );
    
    
    includesstr = {};
    
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinternal%sceres"', slash, slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinternal%sceres%sminiglog"', slash, slash, slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinternal"', slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%sceres-solver%sinclude"', slash, slash)];
    includesstr{end+1} = ['-I"' cpath sprintf('%seigen3"', slash)];
      
    h = msgbox('Ceres is currently compiling and linking. This process can take up to 15 minutes depending on your system. Please be patient...','Please wait...');     
    drawnow;
    
    mex( includesstr{:}, 'CXXOPTIMFLAGS=''-O2''', '-DMAX_LOG_LEVEL=1', sprintf('%sceresd2d.cpp', cpath), ccfileListCeres{:} ,'-outdir', cpath);
    
    close(h);
    
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
    fileList = [fileList; getAllFiles(nextDir)];  %#ok Recursively call getAllFiles
  end

end


