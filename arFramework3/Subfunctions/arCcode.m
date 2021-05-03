% cstr = arCcode(T, matlab_version)
% 
%   T                 any kind of symbolic expression
% 
%   matlab_version    Matlab version e.g. 9.5, use ver('MATLAB')
% 
%   This function generalizes matlab's ccode function and keeps the
%   functionality indpendent of Matlab verions
% 
%   The function was earlier called ccode2 and it was part of arCompileAll
%   It was splitted in order to easier test version dependencies of this
%   conversion step
% 
%   Attention: The resulting might terms seem rather strange, but they are
%   produced in this way by ccode. See arCompileAll how to process them.
% 
% Example:
% T = ar.model(1).condition(1).sym.fv;
% matVer = ver('MATLAB');
% matlab_version = str2double(matVer.Version);
% cstr = arCcode(T, matlab_version)
% cvar = 'fv'
% cstr = strrep(cstr, 't0', [cvar '[0]']);
% cstr = strrep(cstr, '][0]', ']');
% cstr = strrep(cstr, 'T[', [cvar '[']);
% cstr = strrep(cstr, 'udata_splines__', 'data->splines, data->splineIndices');

function cstr = arCcode(T, matlab_version)

% If this matrix or value is empty, do not attempt to generate C-code
if (numel(T) == 0)
    cstr = '';
    return;
end

try
    T = arMyStr2Sym( replaceDerivative( char(T) ) );
catch
    char(T)
    error('Failure attempting to replace derivatives in');
end

% R2015b compatibility fix
if(matlab_version>=8.6)
    sym_str = arSym2str(T);
    if all(strcmp('0',sym_str))
        cstr = char;
    else
        cstr = ccode(T);
    end
else
    cstr = ccode(T);
end

% R2015b compatibility fix
if(matlab_version>=9.8)
    cstr = regexprep(cstr,'t(\d+) =','T[$1][0] =');
end

% Fixes bug with ccode for single line syms
if ( length( cstr ) > 8 )
    if ( strcmp(cstr(1:9), '  _assign') )
        cstr = strrep( [cstr(1:end-2) ';'], '_assign(t0,', 't0=' );
    end
end



function str = replaceDerivative( str )
% Pattern that matches the derivatives D([#], func)(args)
pattern = 'D[\(][\[](\d+)[\]][\,](\s*)(\w*)[\)][\(]';

% Transform D([#], name)(args) => Dname(args, %d)
% We need to explicitly match up the brackets, since regexps can't do
% this reliably, we have to loop over it.
[tokens,startloc,endloc] = regexp( str, pattern, 'tokens' );

% Generate replacement blocks
oldString = cell( 1, numel( startloc ) );
to = cell( 1, numel( startloc ) );
for i = 1 : numel( startloc )
    call = str(startloc(i):endloc(i));
    block = quickScan( str(endloc(i)+1:end) );
    oldString{i} = [ call block ];
    to{i} = sprintf( 'D%s(%s, %s)', tokens{i}{3}, block(1:end-1), tokens{i}{1} );
end

% Perform the replacements
for i = 1 : numel( oldString )
    str = strrep( str, oldString{i}, to{i} );
end

% Pattern which matches the other derivative structure
pattern2 = 'diff[\(]([\w\(\)]+)[\(]([\[\]\(\)\^\/\*\+\-\.\s,\w\d]*)[\)][,](\s)([\[\]\w]*)[\)]';

% Performs regexprep which transforms diff(name(args), args(#)) => Dname(args, #)
[~,~,~,total,matches]=regexp(str, pattern2);

for jm = 1 : numel(matches)
    indep = matches{jm}{4};
    variables = strtrim(strsplit(matches{jm}{2}, ','));
    
    % Find which variable we're deriving w.r.t. to
    ID = find(strcmp(variables, indep));
    
    % Write the new string
    fNew = sprintf( 'D%s(%s, %d)', matches{jm}{1}, matches{jm}{2}, ID );
    
    str = strrep( str, total{jm}, fNew );
end
% because of compatibility with MATLAB2020a and later versions
tmp = ver('MATLAB');
if ~verLessThan('matlab', '9.8')
    str = strrep(str,';',',');
end


function strG = quickScan( str )
    c = 1;
    depth = 0;
    while ( c < numel( str ) )
        if ( str(c) == '(' )
            depth = depth + 1;
        end
        if ( str(c) == ')' )
            depth = depth - 1;
            if ( depth < 0 )
                strG = str(1:c);
                return
            end
        end
        c = c + 1;
    end
    strG = str;
