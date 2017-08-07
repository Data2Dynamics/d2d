% Read a line and parse it. Make sure the entire line is read.
function [C, fid] = arTextScan( fid, varargin )

    tmp = ''; C{1} = {};
    if ( isstruct( fid ) )
        
        if ~isfield( fid, 'nlines' )
            fid.nlines = 0;
        end
        
        % If we are passed a struct, then it is a string with position
        % (parsing from text)
        while ( (fid.pos < numel( fid.str ))&&(isempty(tmp)||isempty(C{1})) )
            % Find next newline
            newlines = regexp(fid.str, '\n|\r\n|\r');
            nl = newlines( find( newlines > fid.pos, 1 ) );
            if ( isempty( nl ) )
                nl = numel( fid.str ) + 1;
            end
            
            % Grab string and advance line
            tmp = fid.str( fid.pos : nl - 1 );
            [C, position] = textscan(tmp, varargin{:});
            
            fid.pos = min( [nl + 1, fid.pos + position] );
            fid.nlines = fid.nlines + 1;
        end
    else
        % If we are passed a file handle, it is business as usual
        while ~feof(fid)&&(isempty(tmp)||isempty(C{1}))
            tmp = fgets(fid);

            C = textscan(tmp, varargin{:} );
        end
    end
    
    if isempty(tmp)
        C = textscan(' ', varargin{:} );
    end
    
end