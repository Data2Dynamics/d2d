% function arValidateInput( C, stage, [varargin] )
% Subroutine used for error checking in the parser
%
% C         indicates a cell array of strings
% stage     should have an identifier for the parsing stage we are in
% varargin  names of the expected function arguments

function arValidateInput( C, stage, varargin )
    for a = 1 : length(varargin)
        if ( iscell(C{a}) )
            try
                C{a} = C{a}{1};
            catch
                if isempty( C{1} )
                    error( 'Did not find data for parsing stage %s. Did you leave a " bracket open?', stage );
                else
                    try
                        str = sprintf( '%s ', C{1} );
                    catch
                        str = '';
                    end
                    
                    error( 'Malformed input during parsing %s at %s', stage, str );
                end
            end
        end 
        if ( isempty( C{a} ) )
            format = getFormat(varargin);
            error( 'Missing %s in %s specification. Correct format is: %s. Supplied was: %s.', varargin{a}, upper(stage), format, fullString(C) );
        end       
        if ( isnan( C{a} ) )
            format = getFormat(varargin);
            error( 'Missing numeric %s parameter in %s specification. Correct format is: %s. Supplied was: %s.', varargin{a}, upper(stage), format, fullString(C) );
        end        
    end
end

function str = fullString(C)
    try
        str = strcat( C{:} );
        if ( isempty( str ) )
            str = '';
        end
    catch
        str = '< Failed to obtain line in which error occurs >';
    end
end

function str = getFormat(varargin)
    str = sprintf( '[%s], ', varargin{1}{:} );
    str = str(1:end-2);
end