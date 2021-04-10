% function arCompile([forceFullCompile, forceCompileLast, debug_mode, source_dir])
% Compile CVODES c-functions
%
%   forceFullCompile:   recompile all objects files     [false]
%   forceCompileLast:   only recompile mex-file         [false]
%   debug_mode:         exclude precompiled objects     [false]
%   source_dir:         external source directory       []
%
% OR
%
% arCompile(ar, forceFullCompile, forceCompileLast, debug_mode, source_dir)
%   ar:                 d2d model/data structure


function arCompile(varargin)

if(nargin==0 || ~isstruct(varargin{1}))
    global ar %#ok<TLEV>
    
    if(isempty(ar))
        error('please initialize by arInit')
    end
else
    ar = varargin{1};
    if(nargin>1)
        varargin = varargin(2:end);
    else
        varargin = {};
    end
end

enableRootfinding = 0;
if ( isfield( ar.config, 'fastEquilibration' ) )
    enableRootfinding = ar.config.fastEquilibration;
end

includeBLAS = 0;
includeLAPACK = 0;
if ( enableRootfinding )
    includeBLAS     = 1;
    includeLAPACK   = 1;
end

global arOutputLevel;
if isempty( arOutputLevel )
    arOutputLevel = 2;
end

if ( arOutputLevel > 2 )
    verbose{1} = '-v';
else
    verbose = {};
end

if ( arOutputLevel > 3 )
    if ( ispc )
        % Full debug output on windows (build a release build, but with
        % debugging symbols. This allows attaching the MSVC debugger)
        % W4        = Warning level 4, produces all relevant warnings
        % Zi        = Put debug info in obj files
        % DEBUG     = Build PDB files for use with the MSVC debugger
        % OPT:REF   = Only include used functions (debug disables this normally)
        % OPT:ICF   = Enable COMDAT folding (debug disables this normally)
        verbose{2} = 'COMPFLAGS="$COMPFLAGS -W4 -Zi"';
        verbose{3} = 'LINKFLAGS="$LINKFLAGS -DEBUG -OPT:REF -OPT:ICF"';
    end
end

mexopt = {'-largeArrayDims'};
if ( includeBLAS )
    mexopt{end+1} = '-lmwblas';
end
if ( includeLAPACK )
    mexopt{end+1} = '-lmwlapack';
end
if ( enableRootfinding )
    mexopt{end+1} = '-DROOT_FINDING';
    ar.config.C_rootfinding = 1;
end
if isfield( ar.config, 'defines' )
    mexopt = union( mexopt, ar.config.defines );
end

arFprintf(1, 'Compiling files...');

usePool = exist('gcp','file')>0 && ~isempty(gcp('nocreate'));

if(~isempty(varargin))
    forceFullCompile = varargin{1};
else
    forceFullCompile = false;
end
if(length(varargin)>1)
    forceCompileLast = varargin{2};
else
    forceCompileLast = false;
end
if(length(varargin)>2)
    debug_mode = varargin{3};
else
    debug_mode = false;
end
if(length(varargin)>3)
    source_dir = varargin{4};
else
    source_dir = pwd;
end

arFprintf(2, '\n');

if(length(which('arClusterCompiledHook.m','-all'))>1)
    warning('arClusterCompiledHook.m is found multiple times which can cause compilation errors. Check your matlab path.');
end

ar_path = fileparts(which('arInit.m'));
% sundials_path = [ar_path filesep '/sundials-2.5.0/' filesep]; % sundials 2.5.0
sundials_path = [ar_path filesep 'ThirdParty' filesep 'sundials-2.6.1' filesep]; % sundials 2.6.1
KLU_path = [ar_path filesep 'ThirdParty' filesep 'KLU-1.2.1'  filesep]; % KLU of suitesparse 4.2.1
compiled_cluster_path = fileparts(which('arClusterCompiledHook.m'));

% compile directory
if(~exist(['./Compiled/' ar.info.c_version_code '/' mexext], 'dir'))
    mkdir(['./Compiled/' ar.info.c_version_code '/' mexext])
end

%% include directories
includesstr = {};

% CVODES
includesstr{end+1} = ['-I"' sundials_path 'include"'];
includesstr{end+1} = ['-I"' sundials_path 'src/cvodes"'];

% KLU
includesstr{end+1} = ['-I"' KLU_path 'SuiteSparse_config"'];
includesstr{end+1} = ['-I"' KLU_path 'KLU/Include"'];
includesstr{end+1} = ['-I"' KLU_path 'KLU/Source"'];
includesstr{end+1} = ['-I"' KLU_path 'AMD/Include"'];
includesstr{end+1} = ['-I"' KLU_path 'AMD/Source"'];
includesstr{end+1} = ['-I"' KLU_path 'BTF/Include"'];
includesstr{end+1} = ['-I"' KLU_path 'BTF/Source"'];
includesstr{end+1} = ['-I"' KLU_path 'COLAMD/Include"'];
includesstr{end+1} = ['-I"' KLU_path 'COLAMD/Source"'];

% arFramework3
includesstr{end+1} = ['-I"' pwd '/Compiled/' ar.info.c_version_code '"'];
includesstr{end+1} = ['-I"' source_dir '/Compiled/' ar.info.c_version_code '"'];
includesstr{end+1} = ['-I"' ar_path '"'];
includesstr{end+1} = ['-I"' ar_path '/Ccode"'];
includesstr{end+1} = ['-I/usr/local/include'];

% Allow interrupt hooks
if ( ar.config.instantaneous_termination )
    includesstr{end+1} = '-lut';
    mexopt{end+1} = '-DALLOW_INTERRUPTS';
end

% Disable code optimization (for code which takes a long time to compile)
if isfield( ar.config, 'no_optimization' )
    if ar.config.no_optimization == 1  || ischar(ar.config.no_optimization)
        
        if ischar(ar.config.no_optimization)  % interpret the config as compiler flag
            mexopt{end+1} = ar.config.no_optimization;
        elseif ispc  % windows
            mexopt{end+1} = 'COMPFLAGS=/Od';
        else
            mexopt{end+1} = 'COMPFLAGS=''-O0''';
        end
        
        fprintf('Optimization disabled\n' );
        warning('ar.config.no_optimization=1: Please carefully check, whether code-optimization really is switched off, e.g. by debugging below where mex is called.')
        fprintf('Default options can overwrite the flag set here (type cc = mex.getCompilerConfigurations). In this case change the mexopt file (see cc.MexOpt). \n');
        fprintf('Microsoft Visual C++ 2013 Professional (C) requires the following flag: %s\n', 'COMPFLAGS=/Od');
    end
end

if(~isempty(compiled_cluster_path))
    includesstr{end+1} = ['-I"' compiled_cluster_path '/' ar.info.c_version_code '"'];
end

c_version_code = ar.info.c_version_code;
objectsstr = {};
%% pre-compile KLU sources

% source files
sourcesKLU = {
    'SuiteSparse_config/SuiteSparse_config.c';
    'KLU/Source/klu.c';
    'KLU/Source/klu_analyze.c';
    'KLU/Source/klu_analyze_given.c';
    'KLU/Source/klu_defaults.c';
    'KLU/Source/klu_diagnostics.c';
    'KLU/Source/klu_dump.c';
    'KLU/Source/klu_extract.c';
    'KLU/Source/klu_factor.c';
    'KLU/Source/klu_free_numeric.c';
    'KLU/Source/klu_free_symbolic.c';
    'KLU/Source/klu_kernel.c';
    'KLU/Source/klu_memory.c';
    'KLU/Source/klu_refactor.c';
    'KLU/Source/klu_scale.c';
    'KLU/Source/klu_solve.c';
    'KLU/Source/klu_sort.c';
    'KLU/Source/klu_tsolve.c';
    'AMD/Source/amd_1.c';
    'AMD/Source/amd_2.c';
    'AMD/Source/amd_aat.c';
    'AMD/Source/amd_control.c';
    'AMD/Source/amd_defaults.c';
    'AMD/Source/amd_dump.c';
    'AMD/Source/amd_global.c';
    'AMD/Source/amd_info.c';
    'AMD/Source/amd_order.c';
    'AMD/Source/amd_post_tree.c';
    'AMD/Source/amd_postorder.c';
    'AMD/Source/amd_preprocess.c';
    'AMD/Source/amd_valid.c';
    'BTF/Source/btf_maxtrans.c';
    'BTF/Source/btf_order.c';
    'BTF/Source/btf_strongcomp.c';
    'COLAMD/Source/colamd.c';
    'COLAMD/Source/colamd_global.c';
    };
sourcesstrKLU = '';
for j=1:length(sourcesKLU)
    sourcesstrKLU = strcat(sourcesstrKLU, [' "' KLU_path sourcesKLU{j} '"']);
end

objectsKLU = {
    'SuiteSparse_config.o';
    'klu.o';
    'klu_analyze.o';
    'klu_analyze_given.o';
    'klu_defaults.o';
    'klu_diagnostics.o';
    'klu_dump.o';
    'klu_extract.o';
    'klu_factor.o';
    'klu_free_numeric.o';
    'klu_free_symbolic.o';
    'klu_kernel.o';
    'klu_memory.o';
    'klu_refactor.o';
    'klu_scale.o';
    'klu_solve.o';
    'klu_sort.o';
    'klu_tsolve.o';
    'amd_1.o';
    'amd_2.o';
    'amd_aat.o';
    'amd_control.o';
    'amd_defaults.o';
    'amd_dump.o';
    'amd_global.o';
    'amd_info.o';
    'amd_order.o';
    'amd_post_tree.o';
    'amd_postorder.o';
    'amd_preprocess.o';
    'amd_valid.o';
    'btf_maxtrans.o';
    'btf_order.o';
    'btf_strongcomp.o';
    'colamd.o';
    'colamd_global.o';
    };
if(ispc)
    objectsKLU = strrep(objectsKLU, '.o', '.obj');
end

for j=1:length(objectsKLU)
    objectsstr = [objectsstr {['./Compiled/' ar.info.c_version_code '/' mexext '/' objectsKLU{j}]}]; %#ok<AGROW>
end


% compile
if(usePool)
    parfor j=1:length(sourcesKLU)
        if(~exist(['Compiled/' c_version_code '/' mexext '/' objectsKLU{j}], 'file') || (forceFullCompile==1))
            mex('-c',verbose{:},mexopt{:}, '-DNTIMER', '-outdir', ['Compiled/' c_version_code '/' mexext '/'], ...
                includesstr{:}, [KLU_path sourcesKLU{j}]); %#ok<PFBNS>
            arFprintf(2, 'compiling KLU(%s)...done\n', objectsKLU{j});
        else
            arFprintf(2, 'compiling KLU(%s)...skipped\n', objectsKLU{j});
        end
    end
else
    for j=1:length(sourcesKLU)
        if(~exist(['Compiled/' c_version_code '/' mexext '/' objectsKLU{j}], 'file') || (forceFullCompile==1))
            mex('-c',verbose{:},mexopt{:}, '-DNTIMER', '-outdir', ['Compiled/' c_version_code '/' mexext '/'], ...
                includesstr{:}, [KLU_path sourcesKLU{j}]);
            arFprintf(2, 'compiling KLU(%s)...done\n', objectsKLU{j});
        else
            arFprintf(2, 'compiling KLU(%s)...skipped\n', objectsKLU{j});
        end
    end
end

%% pre-compile CVODES sources

% source files
sources = {
    'src/cvodes/cvodes_band.c';
    'src/cvodes/cvodes_bandpre.c';
    'src/cvodes/cvodes_bbdpre.c';
    'src/cvodes/cvodes_direct.c';
    'src/cvodes/cvodes_dense.c';
    'src/cvodes/cvodes_klu.c';
    'src/cvodes/cvodes_sparse.c';
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
    'src/sundials/sundials_sparse.c';
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
    'cvodes_klu.o';
    'cvodes_sparse.o';
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
    'sundials_sparse.o';
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
if ( enableRootfinding )
    objects{end+1} = 'inverseC.o';
end
if(ispc)
    objects = strrep(objects, '.o', '.obj');
end

for j=1:length(objects)
    objectsstr = [objectsstr {['./Compiled/' ar.info.c_version_code '/' mexext '/' objects{j}]}]; %#ok<AGROW>
end

% compile
if(usePool)
    parfor j=1:length(sources)
        if(~exist(['Compiled/' c_version_code '/' mexext '/' objects{j}], 'file') || (forceFullCompile==1))
            mex('-c',verbose{:},mexopt{:}, '-outdir', ['Compiled/' c_version_code '/' mexext '/'], ...
                includesstr{:}, [sundials_path sources{j}]); %#ok<PFBNS>
            arFprintf(2, 'compiling CVODES(%s)...done\n', objects{j});
        else
            arFprintf(2, 'compiling CVODES(%s)...skipped\n', objects{j});
        end
    end
else
    for j=1:length(sources)
        if(~exist(['Compiled/' c_version_code '/' mexext '/' objects{j}], 'file') || (forceFullCompile==1))
            mex('-c',verbose{:},mexopt{:}, '-outdir', ['Compiled/' c_version_code '/' mexext '/'], ...
                includesstr{:}, [sundials_path sources{j}]);
            arFprintf(2, 'compiling CVODES(%s)...done\n', objects{j});
        else
            arFprintf(2, 'compiling CVODES(%s)...skipped\n', objects{j});
        end
    end
end

%% pre-compile input functions
if(~ispc)
    objects_inp = ['./Compiled/' ar.info.c_version_code '/' mexext '/arInputFunctionsC.o'];
else
    objects_inp = ['./Compiled/' ar.info.c_version_code '/' mexext '/arInputFunctionsC.obj'];
end

if(~exist(objects_inp, 'file') || forceFullCompile)
    mex('-c',verbose{:},mexopt{:},'-outdir',['Compiled/' ar.info.c_version_code '/' mexext '/'], ...
        includesstr{:}, which('arInputFunctionsC.c'));
    arFprintf(2, 'compiling input functions...done\n');
else
    arFprintf(2, 'compiling input functions...skipped\n');
end

if ( enableRootfinding )
    if(~exist(objects_inp, 'file') || forceFullCompile)
        mex('-c',verbose{:},mexopt{:},'-outdir',['Compiled/' ar.info.c_version_code '/' mexext '/'], ...
            includesstr{:}, which('inverseC.c'));
        arFprintf(2, 'compiling rootfinding functions...done\n');
    else
        arFprintf(2, 'compiling rootfinding functions...skipped\n');
    end
end


% TODO I don't know why this gives a link error ... ?
% if(~debug_mode)
%     objectsstr = [objectsstr {objects_inp}];
% end

%% pre-compile conditions
objects_con = {};
file_con = {};
ms = [];
cs = [];
for jm = 1:length(ar.model)
    for jc = 1:length(ar.model(jm).condition)
        if(~ispc)
            objects_con{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).condition(jc).fkt '.o']; %#ok<AGROW>
        else
            objects_con{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).condition(jc).fkt '.obj']; %#ok<AGROW>
        end
        file_con{end+1} = [ar.model(jm).condition(jc).fkt '.c']; %#ok<AGROW>
        ms(end+1) = jm; %#ok<AGROW>
        cs(end+1) = jc; %#ok<AGROW>
    end
end

if(usePool)
    parfor j=1:length(objects_con)
        if(~exist(objects_con{j}, 'file') || forceFullCompile)
            if(isempty(compiled_cluster_path))
                mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_con{j}]);  %#ok<PFBNS>
            else
                mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_con{j}]);
            end
            arFprintf(2, 'compiling condition m%i c%i, %s...done\n', ms(j), cs(j), file_con{j});
        else
            arFprintf(2, 'compiling condition m%i c%i, %s...skipped\n', ms(j), cs(j), file_con{j});
        end
    end
else
    for j=1:length(objects_con)
        if(~exist(objects_con{j}, 'file') || forceFullCompile)
            if(isempty(compiled_cluster_path))
                mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_con{j}]);
            else
                mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_con{j}]);
            end
            arFprintf(2, 'compiling condition m%i c%i, %s...done\n', ms(j), cs(j), file_con{j});
        else
            arFprintf(2, 'compiling condition m%i c%i, %s...skipped\n', ms(j), cs(j), file_con{j});
        end
    end
end

if(~debug_mode)
    objectsstr = [objectsstr objects_con];
end

%% pre-compile data
if(isfield(ar.model, 'data'))
    objects_dat = {};
    file_dat = {};
    ms = [];
    ds = [];
    
    for jm = 1:length(ar.model)
        for jd = 1:length(ar.model(jm).data)
            if(~ispc)
                objects_dat{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).data(jd).fkt '.o']; %#ok<AGROW>
            else
                objects_dat{end+1} = ['./Compiled/' ar.info.c_version_code '/' mexext '/' ar.model(jm).data(jd).fkt '.obj']; %#ok<AGROW>
            end
            file_dat{end+1} = [ar.model(jm).data(jd).fkt '.c']; %#ok<AGROW>
            ms(end+1) = jm; %#ok<AGROW>
            ds(end+1) = jd; %#ok<AGROW>
        end
    end
    
    if(usePool)
        parfor j=1:length(objects_dat)
            if(~exist(objects_dat{j}, 'file') || forceFullCompile)
                if(isempty(compiled_cluster_path))
                    mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_dat{j}]);  %#ok<PFBNS>
                else
                    mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_dat{j}]);
                end
                arFprintf(2, 'compiling data m%i d%i, %s...done\n', ms(j), ds(j), file_dat{j});
            else
                arFprintf(2, 'compiling data m%i d%i, %s...skipped\n', ms(j), ds(j), file_dat{j});
            end
        end
    else
        for j=1:length(objects_dat)
            if(~exist(objects_dat{j}, 'file') || forceFullCompile)
                if(isempty(compiled_cluster_path))
                    mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_dat{j}]);
                else
                    mex('-c',verbose{:},mexopt{:},'-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_dat{j}]);
                end
                arFprintf(2, 'compiling data m%i d%i, %s...done\n', ms(j), ds(j), file_dat{j});
            else
                arFprintf(2, 'compiling data m%i d%i, %s...skipped\n', ms(j), ds(j), file_dat{j});
            end
        end
    end
    
    if(~debug_mode)
        objectsstr = [objectsstr objects_dat];
    end
end

includesstr=strrep(includesstr,'"', '');
objectsstr = unique( objectsstr, 'stable' );

%% compile and link main mex file
if(~exist([ar.fkt '.' mexext],'file') || forceFullCompile || forceCompileLast)
    if(~ispc)
        % parallel code using POSIX threads for Unix type OS
        chunkSize = getChunkSize( objectsstr );
        chunks = ceil( numel( objectsstr ) / chunkSize );
        if ( chunks > 1 )
            % Link files in chunks to avoid exceeding maximum length (used
            % when there are so many conditions that the command line
            % becomes too big)
            fprintf( 'File string too long for command line. Compiling %d libraries first ...\n', chunks );
            libNames = cell(1, chunks);
            for a = 1 : chunks
                fprintf( 'Assembling chunk %d/%d into library\n', a, chunks );
                libNames{a} = fullfile('Compiled', c_version_code, [ 'lib' num2str(a) '_' ar.fkt(end-31:end) '.o' ] );
                curChunk = cellfun(@(st)[st, ' '], objectsstr( (a-1)*chunkSize+1 : min(a*chunkSize, numel(objectsstr)) ), 'UniformOutput', false);
                system(['ld -r ' [ curChunk{:} ] '-o ' libNames{a}]);
            end
            
            fprintf( 'Compiling system\n' );
            mex(mexopt{:},verbose{:},'-output', ar.fkt, includesstr{:}, '-DHAS_PTHREAD=1', ...
                sprintf('-DNMAXTHREADS=%i', ar.config.nMaxThreads), ...
                which('udata.c'), which('arSimuCalc.c'), libNames{:});
        else
            mex(mexopt{:},verbose{:},'-output', ar.fkt, includesstr{:}, '-DHAS_PTHREAD=1', ...
                sprintf('-DNMAXTHREADS=%i', ar.config.nMaxThreads), ...
                which('udata.c'), which('arSimuCalc.c'), objectsstr{:});
        end
    else
        chunkSize = getChunkSize( objectsstr );
        chunks = ceil( numel( objectsstr ) / chunkSize );
        
        % parallel code using POSIX threads (pthread-win32) for Windows type OS
        includesstr{end+1} = ['-I"' ar_path '\ThirdParty\pthreads-w32_2.9.1\include"'];
        includesstr{end+1} = ['-L"' ar_path '\ThirdParty\pthreads-w32_2.9.1\lib\' mexext '"'];
        includesstr{end+1} = ['-lpthreadVC2'];
        
        if ( chunks > 1 )
            cCfg = mex.getCompilerConfigurations('C');
            libloc = [cCfg.Location, '/VC/bin']; % Add directory containing lib.exe to the path (this could be further generalized)

            curEnv = getenv('path');
            if ( isempty( strfind( curEnv, libloc ) ) )
                setenv( 'path', [curEnv ';' libloc] );
            end

            chunkSize = getChunkSize( objectsstr, 7000 ); % 8192 letter limit in Command and do not want to use powershell -.-
            chunks = ceil( numel( objectsstr ) / chunkSize );
            
            % Link files in chunks to avoid exceeding maximum length (used
            % when there are so many conditions that the command lib becomes too big)
            fprintf( 'File string too long for command line. Compiling %d libraries first ...\n', chunks );
            libNames = cell(1, chunks);

            for a = 1 : chunks
               fprintf( 'Assembling chunk %d/%d into library\n', a, chunks );
               libNames{a} = fullfile('Compiled', c_version_code, [ 'lib' num2str(a) '_' ar.fkt(end-31:end) '.lib' ] );
               curChunk = cellfun(@(st)[st, ' '], objectsstr( (a-1)*chunkSize+1 : min(a*chunkSize, numel(objectsstr)) ), 'UniformOutput', false);               
               system(['vcvars32 & lib /out:' libNames{a} ' ' [ curChunk{:} ] ] );
            end
            fprintf( 'Linking chunks\n' );
            mex(mexopt{:},verbose{:},'-output', ar.fkt, includesstr{:}, '-DHAS_PTHREAD=1', ...
                           sprintf('-DNMAXTHREADS=%i', ar.config.nMaxThreads), ...
                           which('udata.c'), which('arSimuCalc.c'), libNames{:});
        else
            mex(mexopt{:},verbose{:},'-output', ar.fkt, includesstr{:}, '-DHAS_PTHREAD=1', ...
            sprintf('-DNMAXTHREADS=%i', ar.config.nMaxThreads), ...
            which('udata.c'), which('arSimuCalc.c'), objectsstr{:});
        end
    end
    arFprintf(2, 'compiling and linking %s...done\n', ar.fkt);
else
    arFprintf(2, 'compiling and linking %s...skipped\n', ar.fkt);
end

%% refresh file cache
rehash

function chunkSize = getChunkSize( objectsstr, altmax )
if ( nargin < 2 )
    [~, maxArgSize]=system('getconf ARG_MAX');
    maxArgSize = str2num(maxArgSize); %#ok
    if ( isempty( maxArgSize ) )
        maxArgSize = 121071; % Limit is 131071 for windows, but need to have some buffer. %262144;
    end
else
    maxArgSize = altmax;
end

% Take a chunk size with a reasonably safe margin
maxLen = max( cellfun( @numel, objectsstr ) ) + numel(pwd);
chunkSize = floor( (0.5 * maxArgSize) / maxLen );

% Empirically, the maxfile limit seems to be 500 on my system. Not sure
% how to find this number in a platform independent manner.
chunkSize = min( [chunkSize, 500] );
