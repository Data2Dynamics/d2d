% Read a line and parse it. Make sure the entire line is read.
function C = arTextScan( fid, varargin )

    tmp = '';
    while ~feof(fid)&&isempty(tmp)
        tmp = fgets(fid);
        
        if ( min(isspace(tmp)) == 1 )
            tmp = '';
        end
    end

    if ~isempty(tmp)
        C = textscan(tmp, varargin{:} );
    else
        C = textscan(' ', varargin{:} );
    end
    