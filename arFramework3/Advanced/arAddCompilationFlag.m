% arAddCompilationFlag( flag )
%
% adds flag to ar.config.defines
%
% flag - string, name of status
%
% arAddCompilationFlat('1')
% arAddCompilationFlat('Log10')

function arAddCompilationFlag( flag )

	global ar;
    if ~isfield( ar.config, 'defines' )
        ar.config.defines = { flag };
    else
        ar.config.defines = union( ar.config.defines, { flag } );
    end
    fprintf( 'Compilation flags: %s\n', sprintf( '%s ', ar.config.defines{:} ) );