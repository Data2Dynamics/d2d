% better subs

function out = arSubs(in, old, new, matlab_version)

if(nargin<4)
    matVer = ver('MATLAB');
    matlab_version = str2double(matVer.Version);   
end

if(~isnumeric(in) && ~isempty(old) && ~isempty(symvar(in)))
    try
        if(matlab_version>=8.1)
            out = subs(in, old(:), new(:));
        else
            out = subs(in, old(:), new(:), 0);
        end
    catch
        % Failure to substitute, provide some info that might help debug
        % the problem; try them one by one and output those that failed
        s{1} = sprintf( 'Error: Model substitution failure in %s: \n\nThe following substitutions failed:\n', char( in ) );
        for a = 1 : length( old )
            try
                if(matlab_version>=8.1)
                    out = subs(in, old(a), new(a));
                else
                    out = subs(in, old(a), new(a), 0);
                end
            catch ME
                s{end+1} = sprintf( 'Subs [%10s => %5s failed]: %s\n', ...
                    char( old(a) ), char( new(a) ), strtok(ME(1).message, sprintf('\n')) ); %#ok<AGROW>
            end
        end
        s{end+1} = sprintf( '\n\nPlease check substitution errors for clues where the error may be.\n' );
        
        error('%s',s{:});
    end
else
    out = in;
end