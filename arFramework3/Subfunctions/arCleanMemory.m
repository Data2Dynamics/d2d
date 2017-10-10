% Mex files can sometimes leave a huge amount of memory resident.
% This function removes the ones we keep resident, when possible.
function arCleanMemory(clearallcompiled)

    [~, loaded_mexfiles] = inmem('-completenames');
    if ( nargin < 1 )
        clearallcompiled = 0;
    end

    arFprintf( 2, 'Clearing old D2D mex files from memory ' );
    for jm = 1 : numel( loaded_mexfiles )
        [direc, mexName] = fileparts(loaded_mexfiles{jm});

        % Is it one of our files? Then clear it from memory!
        if ( ~isempty( strfind(mexName, 'arSimuCalcFun') ) || (clearallcompiled) )
            try
                arFprintf( 2, '.' );
                oldDir = pwd;
                cd( direc );
                clear( mexName );
                cd( oldDir );
            catch
                warning( 'Could not remove loaded mex file (file locked or directory deleted?)' );
            end
        end
    end
    arFprintf( 2, ' [OK]\n' );
    
end