% Compile CVODES c-functions
%
% arCompile(forceFullCompile)
%
%   forceFullCompile:   recompile all object files      [false]
%
%   This function is based on install_STB written by Radu Serban
%   and installs only the parts of the sundials toolbox.

function arCompile(forceFullCompile, forceCompileLast, debug_mode)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('forceFullCompile','var'))
    forceFullCompile = false;
end
if(~exist('forceCompileLast','var'))
    forceCompileLast = false;
end
if(~exist('debug_mode','var'))
    debug_mode = false;
end

fprintf('\n');

if(~ispc)
    ar_path = strrep(which('arInit.m'),'/arInit.m','');
    sundials_path = [strrep(which('arInit.m'),'/arInit.m','') '/sundials-2.5.0/']; % sundials 2.5.0
else
    ar_path = strrep(which('arInit.m'),'\arInit.m','');
    sundials_path = [strrep(which('arInit.m'),'\arInit.m','') '\sundials-2.5.0\']; % sundials 2.5.0
end

% compile directory
if(~exist([cd '/Compiled/' ar.info.c_version_code '/' mexext], 'dir'))
    mkdir([cd '/Compiled/' ar.info.c_version_code '/' mexext])
end

%% include directories
includesstr = '';

% CVODES
includesstr = strcat(includesstr, [' -I"' sundials_path 'include"']);
includesstr = strcat(includesstr, [' -I"' sundials_path 'src/cvodes"']);

% arFramework3
includesstr = strcat(includesstr, [' -I"' pwd '/Compiled/' ar.info.c_version_code '"']);
includesstr = strcat(includesstr, [' -I"' ar_path '"']);

%% pre-compile CVODES sources
fprintf('compiling CVODES...');

% source files
sources = {
    'src/cvodes/cvodes_band.c';
    'src/cvodes/cvodes_bandpre.c';
    'src/cvodes/cvodes_bbdpre.c';
    'src/cvodes/cvodes_direct.c';
    'src/cvodes/cvodes_dense.c';
    'src/cvodes/cvodes_diag.c';
    'src/cvodes/cvodea.c';
    'src/cvodes/cvodes.c';
    'src/cvodes/cvodes_io.c';
    'src/cvodes/cvodea_io.c';
    'src/cvodes/cvodes_spils.c';
    'src/cvodes/cvodes_spbcgs.c';
    'src/cvodes/cvodes_spgmr.c';
    'src/cvodes/cvodes_sptfqmr.c';
    'src/sundials/sundials_band.c';
    'src/sundials/sundials_dense.c';
    'src/sundials/sundials_iterative.c';
    'src/sundials/sundials_nvector.c';
    'src/sundials/sundials_direct.c';
    'src/sundials/sundials_spbcgs.c';
    'src/sundials/sundials_spgmr.c';
    'src/sundials/sundials_sptfqmr.c';
    'src/sundials/sundials_math.c';
    'src/nvec_ser/nvector_serial.c';
    };
sourcesstr = '';
for j=1:length(sources)
    sourcesstr = strcat(sourcesstr, [' "' sundials_path sources{j} '"']);
end

% objects
objects = {
    'cvodes_band.o';
    'cvodes_bandpre.o';
    'cvodes_bbdpre.o';
    'cvodes_direct.o';
    'cvodes_dense.o';
    'cvodes_diag.o';
    'cvodea.o';
    'cvodes.o';
    'cvodes_io.o';
    'cvodea_io.o';
    'cvodes_spils.o';
    'cvodes_spbcgs.o';
    'cvodes_spgmr.o';
    'cvodes_sptfqmr.o';
    'sundials_band.o';
    'sundials_dense.o';
    'sundials_iterative.o';
    'sundials_nvector.o';
    'sundials_direct.o';
    'sundials_spbcgs.o';
    'sundials_spgmr.o';
    'sundials_sptfqmr.o';
    'sundials_math.o';
    'nvector_serial.o';
    'arInputFunctionsC.o';
    };
if(ispc)
    objects = strrep(objects, '.o', '.obj');
end

objectsstr = '';
for j=1:length(objects)
    objectsstr = strcat(objectsstr, [' Compiled/' ar.info.c_version_code '/' mexext '/' objects{j}]);
end

% compile
for j=1:length(sources)
    if(~exist(['Compiled/' ar.info.c_version_code '/' mexext '/' objects{j}], 'file'))
        fprintf('o');
        eval(['mex -c -outdir Compiled/'  ar.info.c_version_code '/' mexext '/' includesstr ' ' sundials_path sources{j}]);
    end
end
fprintf('...done\n');

%% pre-compile input functions
fprintf('compiling input functions...');

if(~ispc)
    objects_inp = ['./Compiled/' ar.info.c_version_code '/' mexext '/arInputFunctionsC.o'];
else
    objects_inp = ['./Compiled/' ar.info.c_version_code '/' mexext '/arInputFunctionsC.obj'];
end

if(~exist(objects_inp, 'file') || forceFullCompile)
    eval(['mex -c -outdir Compiled/' ar.info.c_version_code '/' mexext '/' includesstr ' ' ar_path '/arInputFunctionsC.c']);
    fprintf('done\n');
else
    fprintf('skipped\n');
end

%% pre-compile conditions
objects_con = {};
for jm = 1:length(ar.model)
    for sc = 1:length(ar.model(jm).condition)
        fprintf('compiling condition m%i c%i, %s...', jm, sc, ar.model(jm).condition(sc).fkt);
        
        if(~ispc)
            objects_con{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).condition(sc).fkt '.o']; %#ok<AGROW>
        else
            objects_con{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).condition(sc).fkt '.obj']; %#ok<AGROW>
        end
        
        if(~exist(objects_con{end}, 'file') || forceFullCompile)
            eval(['mex -c -outdir ./Compiled/' ar.info.c_version_code '/' mexext '/' includesstr ' ' ...
                './Compiled/' ar.info.c_version_code '/' ar.model(jm).condition(sc).fkt '.c']);
            fprintf('done\n');
        else
            fprintf('skipped\n');
        end
    end
end

if(~debug_mode)
    for j=1:length(objects_con)
        objectsstr = strcat(objectsstr, [' ' objects_con{j}]);
    end
end

%% pre-compile data
if(isfield(ar.model, 'data'))
    objects_dat = {};
    
    for jm = 1:length(ar.model)
        for sc = 1:length(ar.model(jm).data)            
            fprintf('compiling data m%i d%i, %s...', jm, sc, ar.model(jm).data(sc).fkt);
            
            if(~ispc)
                objects_dat{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).data(sc).fkt '.o']; %#ok<AGROW>
            else
                objects_dat{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).data(sc).fkt '.obj']; %#ok<AGROW>
            end
            
            if(~exist(objects_dat{end}, 'file') || forceFullCompile)
                eval(['mex -c -outdir ./Compiled/' ar.info.c_version_code '/' mexext '/' includesstr ' ' ...
                    './Compiled/' ar.info.c_version_code '/' ar.model(jm).data(sc).fkt '.c']);
                fprintf('done\n');
            else
                fprintf('skipped\n');
            end
        end
    end
    
    if(~debug_mode)
        for j=1:length(objects_dat)
            objectsstr = strcat(objectsstr, [' ' objects_dat{j}]);
        end
    end
end

%% compile and link main mex file
fprintf('compiling and linking %s...', ar.fkt);
if(~exist([ar.fkt '.' mexext],'file') || forceFullCompile || forceCompileLast)
    outputstr = [' -output ' ar.fkt];
    if(~ispc)
        % parallel code using POSIX threads for Unix type OS
        eval(['mex' outputstr includesstr sprintf(' -DHAS_PTHREAD=1 -DNMAXTHREADS=%i -DHAS_SYSTIME=1 "', ar.config.nMaxThreads) which('arSimuCalc.c') '"' objectsstr]);
    else
        if(exist('C:\Windows\pthreadGC2.dll','file')>0 && exist('C:\Windows\pthreadVC2.dll','file')>0)
            % parallel code using POSIX threads (pthread-win32) for Windows type OS
            includesstr = strcat(includesstr, [' -I"' ar_path '\pthreads-w32_2.9.1\include" -L"' ar_path '\pthreads-w32_2.9.1\lib\' mexext '" -lpthreadVC2']);
            eval(['mex' outputstr includesstr sprintf(' -DHAS_PTHREAD=1 -DNMAXTHREADS=%i "', ar.config.nMaxThreads) which('arSimuCalc.c') '"' objectsstr]);
        else
            % serial code for windows OS
            eval(['mex' outputstr includesstr ' "' which('arSimuCalc.c') '"' objectsstr]);
        end 
    end
    fprintf('done\n');
else
    fprintf('skipped\n');
end

%% refresh file cache
rehash

