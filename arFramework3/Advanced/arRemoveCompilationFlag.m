function arRemoveCompilationFlag( flag )

	global ar;
    if ~isfield( ar.config, 'defines' )
        warning( 'No compilation flags specified' );
    else
        ar.config.defines = setdiff( ar.config.defines, flag );
        fprintf( 'Compilation flags: %s\n', sprintf( '%s ', ar.config.defines{:} ) );
    end
    