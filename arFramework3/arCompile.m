% Compile CVODES c-functions
%
% function arCompile(forceFullCompile, forceCompileLast, debug_mode, source_dir)
%
%   forceFullCompile:   recompile all objects files     [false]
%   forceCompileLast:   only recompile mex-file         [false]
%   debug_mode:         exclude precompiled objects     [false]
%   source_dir:         external source directory       []
% 
% or
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

fprintf('\n');

if(~ispc)
    ar_path = strrep(which('arInit.m'),'/arInit.m','');
    sundials_path = [strrep(which('arInit.m'),'/arInit.m','') '/sundials-2.5.0/']; % sundials 2.5.0
    compiled_cluster_path = strrep(which('arClusterCompiledHook.m'),'/arClusterCompiledHook.m','');
else
    ar_path = strrep(which('arInit.m'),'\arInit.m','');
    sundials_path = [strrep(which('arInit.m'),'\arInit.m','') '\sundials-2.5.0\']; % sundials 2.5.0
    compiled_cluster_path = strrep(which('arClusterCompiledHook.m'),'\arClusterCompiledHook.m','');
end

% compile directory
if(~exist(['./Compiled/' ar.info.c_version_code '/' mexext], 'dir'))
    mkdir(['./Compiled/' ar.info.c_version_code '/' mexext])
end

%% include directories
includesstr = {};

% CVODES
includesstr{end+1} = ['-I"' sundials_path 'include"'];
includesstr{end+1} = ['-I"' sundials_path 'src/cvodes"'];

% arFramework3
includesstr{end+1} = ['-I"' pwd '/Compiled/' ar.info.c_version_code '"'];
includesstr{end+1} = ['-I"' source_dir '/Compiled/' ar.info.c_version_code '"'];
includesstr{end+1} = ['-I"' ar_path '"'];
if(~isempty(compiled_cluster_path))
    includesstr{end+1} = ['-I"' compiled_cluster_path '/' ar.info.c_version_code '"'];
end

c_version_code = ar.info.c_version_code;

%% pre-compile CVODES sources

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

objectsstr = {};
for j=1:length(objects)
    objectsstr = [objectsstr {['./Compiled/' ar.info.c_version_code '/' mexext '/' objects{j}]}]; %#ok<AGROW>
end

% compile
if(usePool)
    parfor j=1:length(sources)
        if(~exist(['Compiled/' c_version_code '/' mexext '/' objects{j}], 'file'))
            mex('-c', '-outdir', ['Compiled/' c_version_code '/' mexext '/'], ...
                includesstr{:}, [sundials_path sources{j}]); %#ok<PFBNS>
            fprintf('compiling CVODES(%s)...done\n', objects{j});
        else
            fprintf('compiling CVODES(%s)...skipped\n', objects{j});
        end
    end
else
    for j=1:length(sources)
        if(~exist(['Compiled/' c_version_code '/' mexext '/' objects{j}], 'file'))
            mex('-c', '-outdir', ['Compiled/' c_version_code '/' mexext '/'], ...
                includesstr{:}, [sundials_path sources{j}]);
            fprintf('compiling CVODES(%s)...done\n', objects{j});
        else
            fprintf('compiling CVODES(%s)...skipped\n', objects{j});
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
    mex('-c','-outdir',['Compiled/' ar.info.c_version_code '/' mexext '/'], ...
        includesstr{:}, [ar_path '/arInputFunctionsC.c']);
    fprintf('compiling input functions...done\n');
else
    fprintf('compiling input functions...skipped\n');
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
                mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_con{j}]);  %#ok<PFBNS>
            else
                mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_con{j}]);
            end
            fprintf('compiling condition m%i c%i, %s...done\n', ms(j), cs(j), file_con{j});
        else
            fprintf('compiling condition m%i c%i, %s...skipped\n', ms(j), cs(j), file_con{j});
        end
    end
else
    for j=1:length(objects_con)
        if(~exist(objects_con{j}, 'file') || forceFullCompile)
            if(isempty(compiled_cluster_path))
                mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_con{j}]);
            else
                mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                    includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_con{j}]);
            end
            fprintf('compiling condition m%i c%i, %s...done\n', ms(j), cs(j), file_con{j});
        else
            fprintf('compiling condition m%i c%i, %s...skipped\n', ms(j), cs(j), file_con{j});
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
                    mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_dat{j}]);  %#ok<PFBNS>
                else
                    mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_dat{j}]);
                end
                fprintf('compiling data m%i d%i, %s...done\n', ms(j), ds(j), file_dat{j});
            else
                fprintf('compiling data m%i d%i, %s...skipped\n', ms(j), ds(j), file_dat{j});
            end
        end
    else
        for j=1:length(objects_dat)
            if(~exist(objects_dat{j}, 'file') || forceFullCompile)
                if(isempty(compiled_cluster_path))
                    mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [source_dir '/Compiled/' c_version_code '/' file_dat{j}]);
                else
                    mex('-c','-outdir',['./Compiled/' c_version_code '/' mexext '/'], ...
                        includesstr{:}, [compiled_cluster_path '/' c_version_code '/' file_dat{j}]);
                end
                fprintf('compiling data m%i d%i, %s...done\n', ms(j), ds(j), file_dat{j});
            else
                fprintf('compiling data m%i d%i, %s...skipped\n', ms(j), ds(j), file_dat{j});
            end
        end
    end
    
    if(~debug_mode)
        objectsstr = [objectsstr objects_dat];
    end
end

includesstr=strrep(includesstr,'"', '');
%% compile and link main mex file
if(~exist([ar.fkt '.' mexext],'file') || forceFullCompile || forceCompileLast)
    if(~ispc)
        % parallel code using POSIX threads for Unix type OS

        mex('-output', ar.fkt, includesstr{:}, '-DHAS_PTHREAD=1', ...
            sprintf('-DNMAXTHREADS=%i', ar.config.nMaxThreads), ...
            which('arSimuCalc.c'), objectsstr{:});
    else
        if(exist('C:\Windows\pthreadGC2.dll','file')>0 && exist('C:\Windows\pthreadVC2.dll','file')>0)
            % parallel code using POSIX threads (pthread-win32) for Windows type OS
            
            includesstr{end+1} = ['-I"' ar_path '\pthreads-w32_2.9.1\include"'];
            includesstr{end+1} = ['-L"' ar_path '\pthreads-w32_2.9.1\lib\' mexext '"'];
            includesstr{end+1} = '-lpthreadVC2';
            
            mex('-output', ar.fkt, includesstr{:}, '-DHAS_PTHREAD=1', ...
                sprintf('-DNMAXTHREADS=%i', ar.config.nMaxThreads), ...
                which('arSimuCalc.c'), objectsstr{:});
        else
            % serial code for windows OS
            
            mex('-output', ar.fkt, includesstr{:}, ...
                which('arSimuCalc.c'), objectsstr{:});
        end 
    end
    fprintf('compiling and linking %s...done\n', ar.fkt);
else
    fprintf('compiling and linking %s...skipped\n', ar.fkt);
end

%% refresh file cache
rehash

