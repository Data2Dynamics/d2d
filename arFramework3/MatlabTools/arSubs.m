% out = arSubs(in, old, new, matlab_version)
% Better Symbolic substitution than subs.
% 
% in    symbolic expression
% old   sym/str variable to be replaced
% new   sym/num variable with which to replace
%
% Example:
%       subs(cos(a)+sin(b),{a,b},{sym('alpha'),2}) returns
%       cos(alpha)+sin(2)
%
% See also SUBS, SYM

function out = arSubs(in, old, new, matlab_version)
if ischar(old)
    error('Variable ''old'' should be symbolic.');
end
if ischar(new)
    error('Variable ''new'' should be symbolic.');
end

if(nargin<4)
    matVer = ver('MATLAB');
    matlab_version = str2double(matVer.Version);   
end

if(~isnumeric(in) && ~isempty(old) && ~isempty(symvar(in)))
    try
        if(matlab_version>=8.1)
            if size(old,1) ~= size(new,1)
                new = new.';
            end
            new = TransformNew(new);
            %subset = logical(old~=new);
            out = subs(in(1), old, new);
        else
            new = TransformNew(new);
            out = subs(in, old(:), new(:), 0);
        end
    catch ME
        % Failure to substitute, provide some info that might help debug
        % the problem; try them one by one and output those that failed
        s{1} = sprintf( 'Error: Model substitution failure in %s (%s): \n\nThe following substitutions failed:\n', char( in ), ME.message );
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
end
function newOut = TransformNew(newIn)
    for i = 1:length(newIn)
       if newIn{i}(1) == '(' && newIn{i}(end) == ')'
           newIn{i}(1) = '';
           newIn{i}(end) = '';
       end
    end
    newOut = newIn;
end