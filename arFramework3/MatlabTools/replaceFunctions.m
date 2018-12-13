% Replace a function call using an output mask.

function [str, arg] = replaceFunctions(str, funcTypes, checkValidity)

    if (nargin < 3)
        checkValidity = 0;
    end

    arg = [];
    str  = char(str);
    stro = str; replaced = 0; %#ok
    for a = 1 : length( funcTypes )
        funcs = findFunc( str, funcTypes{a}{1} );
        argLayout = funcTypes{a}{3};
        
        for b = 1 : length( funcs )
            if ( length( funcs(b).args ) ~= max(argLayout) )
                msg = { 'Invalid number of function argument for function "', ...
                        funcTypes{a}{1}, '" in model or data definition file. Expected ', ...
                        num2str(max(argLayout)), ' got ', num2str( length( funcs(b).args ) ), ...
                        '. Valid syntax would be ', funcTypes{a}{4} };
                error( sprintf( '%s', msg{:} ) ); %#ok
            else
                % Determine what the function should be replaced with;
                % feed the appropriate function arguments and replace it
                % Also making sure to use extra brackets for safety (e.g.
                % 5*(a+b) != 5*a+b)
                try
                    to = sprintf( ['(' funcTypes{a}{2} ')'], funcs(b).args{funcTypes{a}{3}} );
                    str = strrep( str, funcs(b).func, to );
                    replaced = replaced + 1;
                catch
                    msg = { 'Failed to replace function ', funcTypes{a}{1}, ...
                        ' in:', funcs(b).func, 'Please expression check for error.' };
                    error( sprintf( '%s\n', msg{:} ) ); %#ok
                end
            end
        end
        if ( nargout > 1 )
            arg.(funcTypes{a}{1}) = funcs;
        end
    end
    
    % Determine whether we got a valid symbolic expression and optionally
    % simplify it
    try
        if (checkValidity)
            str = char( arMyStr2Sym( str ) );
        end
        % Enable for input function debug purposes
        % if ( replaced > 0 )
        %     disp(sprintf( '%s =>\n\t\t%s', stro, str ));
        % end
    catch
        msg = { 'Failed to obtain valid expression from: ', ...
                str, 'Please check expression for error.' };
        error(sprintf('%s\n', msg{:})); %#ok
        
    end
    
% Replace variables in the arguments of specific functions
function str = replaceWithinFunc(str, func, from, to) %#ok
    funcs = findFunc( str, func );
    for a = 1 : length( funcs )
        str = strrep(str, funcs(a).func, strrep(funcs(a).func,from,to));
    end
    
% Function to scan for specific function name and extract its arguments
function [f] = findFunc( st, funcName )
    loc     = strfind( st, [funcName '('] );
    if ( length(loc) > 0 ) %#ok
        for a = 1 : length( loc )
            f(a) = fetchArgs( st(loc(a):end) ); %#ok
            f(a).fin = f(a).fin + loc(a)-1; %#ok
        end
    else
        f = [];
    end
 
% Function to fetch function arguments
function f = fetchArgs( st )
    commas  	= [];
    cur         = 0;
    brackets    = 0;
    while( brackets == 0 )
        cur = cur + 1;
        if ( cur > length( st ) )
            error( sprintf( 'Malformed input string for argument fetcher: \n%s', st ) ); %#ok
        end
        if ( brackets < 0 )
            error( sprintf( 'Malformed input string for argument fetcher: \n%s', st ) ); %#ok
        end
        if ( st( cur ) == '(' )
            brackets = brackets + 1;
        end
        if ( st( cur ) == ')' )
            brackets = brackets - 1;
        end        
    end
    if ( brackets < 0 )
        error( sprintf( 'Malformed input string for argument fetcher: \n%s', st ) ); %#ok
    end    
    
    f.name = strtrim( st(1:cur-1) );
    stPos = cur;
    
    while( brackets > 0 )
        cur = cur + 1;
        if ( cur > length( st ) )
            error( sprintf( 'Malformed input string for argument fetcher: \n%s', st ) ); %#ok
        end            
        if ( st( cur ) == '(' )
            brackets = brackets + 1;
        end
        if ( st( cur ) == ')' )
            brackets = brackets - 1;
        end
        if ( ( st( cur ) == ',' ) && ( brackets == 1 ) )
            commas(end+1) = cur; %#ok
        end
    end
    
    f.fin    = cur;
        
    list = [stPos, commas, f.fin];
    for b = 1 : length( list ) - 1
        f.args{b} = strtrim( st(list(b)+1:list(b+1)-1) );
    end
    
    f.func = st(1:cur);