% Compile CVODES c-functions
%
% arCompile(forceFullCompile)
%   forceFullCompile:   recompile all object files      [false]
%
%   This function is based on install_STB written by Radu Serban
%   and installs only the parts of the sundials toolbox.

function arCompile(forceFullCompile)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('forceFullCompile','var'))
    forceFullCompile = false;
end

fprintf('\ncompiling...');

if(exist([ar.fkt '.' mexext],'file') && ~forceFullCompile)
    fprintf('skipped\n');
    return
end

ar_path = strrep(which('arInit.m'),'/arInit.m','');

% sundials_path = [strrep(which('arInit.m'),'/arInit.m','') '/sundials-2.4.0/']; % sundials 2.4.0
sundials_path = [strrep(which('arInit.m'),'/arInit.m','') '/sundials-2.5.0/']; % sundials 2.5.0

% include directories
includes = {'include', 'src/cvodes'};
includesstr = '';
for j=1:length(includes)
    includesstr = strcat(includesstr, [' -I"' sundials_path includes{j} '"']);
end
includesstr = strcat(includesstr, [' -I"' pwd '/Compiled"']);
includesstr = strcat(includesstr, [' -I"' ar_path '"']);

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

objectsstr = '';
for j=1:length(objects)
    objectsstr = strcat(objectsstr, [' Compiled/' mexext '/' objects{j}]);
end

% compile

if(~exist([cd '/Compiled/' mexext], 'dir'))
    mkdir([cd '/Compiled/' mexext])
end

% serial code
% outputstr = ' -output arSimuCalcSerial';

% parallel code using POSIX threads
outputstr = [' -output ' ar.fkt];

% pre-compile CVODES sources
for j=1:length(sources)
    if(~exist(['Compiled/' mexext '/' objects{j}], 'file'))
        fprintf('o');
        eval(['mex -c -outdir Compiled/'  mexext '/' includesstr ' ' sundials_path sources{j}]);
    end
end

% pre-compile input functions
eval(['mex -c -outdir Compiled/'  mexext '/' includesstr ' ' ar_path '/arInputFunctionsC.c']);

% pre-compile conditions
sources_con = {};
objects_con = {};

for jm = 1:length(ar.model)
    for sc = 1:length(ar.model(jm).condition)
        sources_con{end+1} = ['./Compiled/' ar.model(jm).condition(sc).fkt '.c']; %#ok<AGROW>
        objects_con{end+1} = ['./Compiled/' mexext '/' ar.model(jm).condition(sc).fkt '.o']; %#ok<AGROW>
    end
end

for j=1:length(objects_con)
    objectsstr = strcat(objectsstr, [' ' objects_con{j}]);
end

for j=1:length(sources_con)
    if(~exist(objects_con{j}, 'file'))
        fprintf('o');
        eval(['mex -c -outdir ./Compiled/' mexext '/' includesstr ' ' sources_con{j}]);
    end
end

% pre-compile data
if(isfield(ar.model, 'data'))
    sources_dat = {};
    objects_dat = {};
    
    for jm = 1:length(ar.model)
        for sc = 1:length(ar.model(jm).data)
            sources_dat{end+1} = ['./Compiled/' ar.model(jm).data(sc).fkt '.c']; %#ok<AGROW>
            objects_dat{end+1} = ['./Compiled/' mexext '/' ar.model(jm).data(sc).fkt '.o']; %#ok<AGROW>
        end
    end
    
    for j=1:length(objects_dat)
        objectsstr = strcat(objectsstr, [' ' objects_dat{j}]);
    end
    
    for j=1:length(sources_dat)
        if(~exist(objects_dat{j}, 'file'))
            fprintf('o');
            eval(['mex -c -outdir ./Compiled/' mexext '/' includesstr ' ' sources_dat{j}]);
        end
    end
end

% compile and link main mex file
fprintf('o');
eval(['mex' outputstr includesstr ' "' which('arSimuCalc.c') '"' objectsstr]);

fprintf('... done\n');

% refresh file cache
rehash

