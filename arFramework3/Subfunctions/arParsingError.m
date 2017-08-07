% arError
function arParsingError( fid, varargin )

    if isstruct(fid)
        [~,f,e]=fileparts(fid.fn);
        fid.nlines
        varargin{1} = [ varargin{1} sprintf( ' in <a href="matlab: opentoline(''%s'', %d)">%s</a>', strrep(fid.fn, '\', '\\'), fid.nlines, [f e] ) ];
        error( varargin{:} );
    else
        error( varargin{:} );
    end
