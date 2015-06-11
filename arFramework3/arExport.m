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
    
    % Grab all filenames in the data folder, so we have something to match against
    fileList = listFiles( {}, 'Data' );
    fileListU = strrep(fileList,'/','_');    
    
    setupFile = sprintf('arInit;\n\n');
    h = waitbar(0);
    for m = 1 : length( models )
        data = {};
        
        % Grab the models
        modelName = [ 'Models/' ar.model(m).name '.def' ];
        modelDescription = sprintf('%s ', ar.model(m).description{:});
        setupFile = sprintf( '\n%s\n%%%s\n\narLoadModel(''%s'');\n', setupFile, modelDescription, ar.model(m).name );
        mName = strrep(ar.model(m).name, '_', '\_');
        if ( ~exist( modelName, 'file' ) )
            error( 'Cannot find model file %s', modelName );
        else
            dir = fetchDir( modelName );
            if ~isempty( dir )
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
                data = union( data, ['Data_' ar.model(m).data(d).name] );
            end
        end
        
        % Do the data copies
        for d = 1 : length( data )
            % Check which files exist and match up (since the data name mapping is non-unique)
            waitbar( d/length(data), h, sprintf( 'Model %s: Exporting data file %d/%d', mName, d, length(data) ) );
            [files, dir] = findMatchingFiles( data{d}, fileList, fileListU );
            
            if ( ~isempty( dir ) )
                warning('off', 'MATLAB:MKDIR:DirectoryExists');
                mkdir( [ directory '/' dir ] );
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
                        setupFile = sprintf( '%sarLoadData(''%s'', ''csv'');\n', setupFile, files{f}(6:end-4) );
                    end
                end
            end
        end 
    end
    
    arCompress;
    save( [ directory '/Results/arStruct.mat' ], 'ar' );
    setupFile = sprintf( '%s\n\narCompileAll;\n\nload(''Results/arStruct.mat'');\narChi2(false);\narPlotY;', setupFile );
    
    fid = fopen( [directory '/Setup.m' ], 'w+' );
    fprintf( fid, '%s', setupFile );
    fclose(fid);
    
    close(h);
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
