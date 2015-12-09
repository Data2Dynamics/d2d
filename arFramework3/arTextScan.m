% Read a line and parse it. Make sure the entire line is read.
function C = arTextScan( fid, varargin )

    tmp = ''; C{1} = {};
    while ~feof(fid)&&(isempty(tmp)||isempty(C{1}))
        tmp = fgets(fid);
        
        C = textscan(tmp, varargin{:} );
    end

    if ~isempty(tmp)
        C = textscan(tmp, varargin{:} );
    else
        C = textscan(' ', varargin{:} );
    end
    