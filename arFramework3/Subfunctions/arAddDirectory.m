% Adds a directory to the directories logged when using arExport.
% This is useful for keeping track of dependencies.

function arAddDirectory( path )

    global ar;    

    if ( nargin < 1 )
        error( 'Too few arguments in call' );
    end
    
    directories = {};
    if ( isfield( ar, 'directories' ) )
        directories = ar.directories;
    end
    
    A = dir( ['*' path] );
    if isempty(A)
        error( 'Directory does not exist' );
    end
    if ( A.isdir == 0 )
        error( 'Path is not a directory' );
    end
    
    ar.directories = union( directories, path );

end