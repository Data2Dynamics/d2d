% Function used to package model file(s) and imported data for distribution.
%
% Traverses the ar structure to see which data sets and model files are  
% actually used and packs those into the directory specified by the user.
% This includes the data files (xls, xlsx, csv), all the def files (model
% and data), a file with the compressed ar structure (contains parameter set 
% and event information) and a little script that compiles the model and 
% loads the ar struct with the parameters and events.
%
% Function arExport(directory, (models))
%
% Directory     target directory for the package
% Models        model indices (optional)
%
% Example:
%    arExport( 'potato' );
%           puts all the files into the subdirectory potato.

function arExport(directory, models)

    global ar;
    
    if nargin < 2
        models = 1 : length( ar.model );
    end
    
    if(~exist('Models','dir'))
        error('folder Models/ does not exist')
    end
   
    if ( isdir( [ pwd '/' directory ] ) )
        fprintf( '*** WARNING: TARGET DIRECTORY %s ALREADY EXISTS\n', directory );
        fprintf( '*** CONTENTS WILL BE ERASED IF YOU PROCEED! ***\n' );
        fprintf( 'Hit any key to continue or CTRL+C to abort\n' );
        pause;
        rmdir([ pwd '/' directory ],'s');
    end
    
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir( directory );
    mkdir( [ directory, '/Models' ] );
    mkdir( [ directory, '/Data' ] );
    mkdir( [ directory, '/Results' ] );
    warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
    if(~isfield(ar.config,'useFitErrorMatrix'))
        ar.config.useFitErrorMatrix = false;
    end
    
    % Grab all filenames in the data folder, so we have something to match against
    fileList = listFiles( {}, 'Data' );
    fileListU = strrep(fileList,'/','_');
    fileListU = strrep(fileListU,'-','_');
    
    subPaths = copyDependencies( directory );
    
    % If mexfiles exist, copy them!
    mx = dir([ar.fkt '.*']); %#ok
    for a = 1 : length( mx );
        copyfile( mx(a).name, directory );
    end
    
    setupFile = sprintf('%% D2D Model Package\n%%   Model author: %s.\n%%   More information: http://data2dynamics.org\n\n%% Initialize D2D framework\n%% Add extra paths\n%s\narInit;\n', ar.config.username, subPaths );
    h = waitbar(0);
    for m = 1 : length( models )
        data = {};
        
        % Grab the models
        modelName = [ 'Models/' ar.model(m).name '.def' ];
        modelDescription = sprintf('%s ', ar.model(m).description{:});
        setupFile = sprintf( '%s\n%% Load model and data for %s.def\n%% %s\narLoadModel(''%s'');\n', setupFile, ar.model(m).name, modelDescription, ar.model(m).name );
        mName = strrep(ar.model(m).name, '_', '\_');
        if ( ~exist( modelName, 'file' ) )
            error( 'Cannot find model file %s', modelName );
        else
            dirn = fetchDir( modelName );
            if ~isempty( dirn )
                warning('off', 'MATLAB:MKDIR:DirectoryExists');
                mkdir( [ directory '/' fetchDir( modelName ) ] );
                warning('on', 'MATLAB:MKDIR:DirectoryExists');            
            end
            copyfile( modelName, [directory '/' modelName] );
        end
        
        % Do we have data?
        if isfield( ar.model(m), 'data' )
            if(~exist('Data','dir'))
                error('folder Data/ does not exist');
            end
            
            % Fetch the names of the files
            for d = 1 : length( ar.model(m).data )
                name = ar.model(m).data(d).name;
                
                % Do we have RANDOM name substitutions? ==> Remove these.
                if ( isfield( ar.model(m).data(d), 'fprand' ) )
                    for jr = 1 : length( ar.model(m).data(d).prand )
                        name = strrep( name, sprintf('_%s%s', ar.model(m).data(d).prand{jr}, ar.model(m).data(d).fprand(jr)), '');
                    end
                end
                data = union( data, ['Data_' name] );
            end
        end
        
        % Do the data copies
        for d = 1 : length( data )
            % Check which files exist and match up (since the data name mapping is non-unique)
            waitbar( d/length(data), h, sprintf( 'Model %s: Exporting data file %d/%d', mName, d, length(data) ) );
            [files, dirn] = findMatchingFiles( data{d}, fileList, fileListU );
            
            if ( ~isempty( dirn ) )
                warning('off', 'MATLAB:MKDIR:DirectoryExists');
                mkdir( [ directory '/' dirn ] );
                warning('on', 'MATLAB:MKDIR:DirectoryExists');
                for f = 1 : length( files )
                    copyfile( files{f}, [directory '/' files{f}] );
                    
                    ff = fliplr(files{f});
                    if ( strfind( ff, fliplr( 'xlsx' ) ) ) 
                        setupFile = sprintf( '%sarLoadData(''%s'');\n', setupFile, files{f}(6:end-5) );
                    end
                    if ( strfind( ff, fliplr( 'xls' ) ) )
                        setupFile = sprintf( '%sarLoadData(''%s'');\n', setupFile, files{f}(6:end-4) );
                    end                    
                    if ( strfind( ff, fliplr( 'csv' ) ) == 1 )
                        setupFile = sprintf( '%sarLoadData(''%s'', [], ''csv'');\n', setupFile, files{f}(6:end-4) );
                    end
                end
            end
        end 
    end
    
    arSaveParOnly( ar, sprintf( '%s/Results/%s', directory, directory ) );
    
    if ( isfield( ar, 'eventLog' ) && (length(ar.eventLog)>0) ) %#ok
        eventCommands = sprintf('%s;\n', ar.eventLog{:} );
        eventCommands = sprintf('%% Setup events\n%s', eventCommands );
    else
        eventCommands = '';
    end
    
    activeDataCommands = activeDatasets();
    
    configFields = { 'rtol', 'atol', 'eq_tol', 'eq_step_factor', 'init_eq_step', 'max_eq_steps', 'maxsteps', 'maxstepsize', 'nFinePoints', 'atolV', 'atolV_Sens', 'ssa_min_tau', 'ssa_runs', 'steady_state_constraint', 'useEvents', 'useFitErrorCorrection', 'useFitErrorMatrix'};
    if(ar.config.useFitErrorMatrix==0)
        configFields = [configFields, 'fiterrors', 'ploterrors'];
    else
        configFields = [configFields, 'fiterrors_matrix', 'ploterrors_matrix'];
    end
        
    configFile      = sprintf( '%% D2D Configuration file\n%s', fieldcode( ar.config, configFields{:} ) );
    
    setupFile = sprintf( '%s\n%% Compile model\narCompileAll;\n\n%% Load configuration\nmodelConfig\n\n%% Load parameters\narLoadPars(''%s'');\n\n%s\n%% Do not fit specific conditions\n%s\n%% Simulate and plot\narSimu(false,true,true);\narChi2(false);\narPlotY;', setupFile, directory, eventCommands, activeDataCommands );
    
    fid = fopen( [directory '/Setup.m' ], 'w+' );
    fprintf( fid, '%s', setupFile );
    fclose(fid);
    
    fid = fopen( [directory '/modelConfig.m' ], 'w+' );
    fprintf( fid, '%s', configFile );
    fclose(fid);    
    
    close(h);
end

function str = activeDatasets()
    global ar;
    
    lines = {};
    for m = 1 : length( ar.model )
        for d = 1 : length( ar.model(m).data )
            name = ar.model(m).data(d).name;
            if ( sum( ar.model(m).data(d).qFit ) == 0 )
                lines{end+1} = sprintf( 'arDisableData(arFindData(''%s'', ''exact'', ''names''));', name ); %#ok<AGROW>
            else
                if ( sum( ar.model(m).data(d).qFit ) < length( ar.model(m).data(d).qFit ) )
                    lines{end+1} = sprintf( 'ar.model(%d).data(arFindData(''%s'', ''exact'')).qFit = [%s];', m, name, sprintf('%d ', ar.model(m).data(d).qFit ) ); %#ok<AGROW>
                end
            end
        end
    end
    str = sprintf( '%s\n', lines{:} );
end

function str = copyDependencies(savepath)
    global ar;

    lines = {};
    if isfield( ar, 'directories' )
        for a = 1 : length( ar.directories )
            warning('off', 'MATLAB:MKDIR:DirectoryExists');
            copyfile( ar.directories{a}, [savepath '/' ar.directories{a}]);
            warning('on', 'MATLAB:MKDIR:DirectoryExists');
            lines{end+1} = sprintf( 'addpath %s;\n', ar.directories{a} ); %#ok<AGROW>
        end
    end
    str = sprintf( '%s', lines{:} );
end

function str = fieldcode( struct, varargin ) %#ok
    global ar;
    
    maxLen = max(cellfun(@length, varargin));
    
    str = [];
    for a = 1 : length( varargin )
        str = sprintf( '%sar.config.%s = %s;\n', str, arExtendStr( varargin{a}, maxLen ), mat2str(ar.config.(varargin{a})) );
    end
end

% save only parameters
function arSaveParOnly(ar2, savepath)

    ar = struct([]);
    ar(1).pLabel = ar2.pLabel;
    ar.p = ar2.p;
    ar.qLog10 = ar2.qLog10;
    ar.qFit = ar2.qFit;
    ar.lb = ar2.lb;
    ar.ub = ar2.ub;
    ar.type = ar2.type;
    ar.mean = ar2.mean;
    ar.std = ar2.std;
    try %#ok<TRYNC>
        ar.chi2fit = ar2.chi2fit;
        ar.ndata = ar2.ndata;
        ar.nprior = ar2.nprior;
    end
    try %#ok<TRYNC>
        ar.ps = ar2.ps;
        ar.ps_errors = ar2.ps_errors;
        ar.chi2s = ar2.chi2s;
        ar.chi2sconstr = ar2.chi2sconstr;
        ar.timing = ar2.timing;
        ar.exitflag = ar2.exitflag;
        ar.fun_evals = ar2.fun_evals;
        ar.ps_start = ar2.ps_start;
        ar.chi2s_start = ar2.chi2s_start;
        ar.chi2sconstr_start = ar2.chi2sconstr_start;
        ar.optim_crit = ar2.optim_crit;

        ar.chi2s_sorted = ar2.chi2s_sorted;
        ar.chi2sconstr_sorted = ar2.chi2sconstr_sorted;
        ar.ps_sorted = ar2.ps_sorted;
        ar.chi2s_start_sorted = ar2.chi2s_start_sorted;
        ar.chi2sconstr_start_sorted = ar2.chi2sconstr_start_sorted;
        ar.ps_start_sorted = ar2.ps_start_sorted;
    end
    if(ar2.config.useFitErrorMatrix == 0)
        ar.config.fiterrors = ar2.config.fiterrors; %#ok<STRNU>
    else
        ar.config.useFitErrorMatrix = true;
        ar.config.fiterrors_matrix = ar2.config.fiterrors_matrix; 
        ar.config.ploterrors_matrix = ar2.config.ploterrors_matrix;
    end
    
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir( savepath );
	warning('on', 'MATLAB:MKDIR:DirectoryExists');    
    save([savepath '/workspace_pars_only.mat'],'ar','-v7.3');
end

function fileList = listFiles( fileList, directory )
    df = dir( directory );
    for f = 1 : length(df)
        if ( df(f).isdir )
            if ( ~strcmp( df(f).name, '.' ) ) && ( ~strcmp( df(f).name, '..' ) )
                fileList = listFiles( fileList, [ directory '/' df(f).name ] );
            end
        else
            fileList = [ fileList ; [ directory '/' df(f).name ] ];%#ok<AGROW>
        end
    end
end

function [matchedFiles, matchedDir, dataFiles] = findMatchingFiles( str, fileList, fileListU )
    id = ismember( fileListU, [str '.def'] );
    matchedFiles = fileList( id );
    matchedDir = {};
    if ( length( matchedFiles ) > 1 )
        warning( 'Found two datasets matching %s. Copying both.', str );
    end
    if ( isempty( matchedFiles ) )
        warning( 'Could not locate data file for %s.', str );
    else
        % Find corresponding xls and csv files
        dataFiles = {};
        dataFiles = union( dataFiles, fileList( ismember( fileListU, [str '.xls'] ) ) );
        dataFiles = union( dataFiles, fileList( ismember( fileListU, [str '.xlsx'] ) ) );
        dataFiles = union( dataFiles, fileList( ismember( fileListU, [str '.csv'] ) ) );
        matchedFiles = union( matchedFiles, dataFiles );
        
        matchedDir   = fetchDir( matchedFiles{1} );
    end
end

function dir = fetchDir( str )
	id = strfind(str,'/');
    if ( ~isempty( id ) )
        dir = str(1:id(end));
    else
        dir = '';
    end
end

