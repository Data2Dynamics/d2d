% Subroutine used for error checking in the parser
% function arValidateInput( C, stage, varargin )
%
% C         indicates a cell array of strings
% stage     should have an identifier for the parsing stage we are in
% varargin  names of the expected function arguments
function arValidateInput( C, stage, varargin )
    for a = 1 : length(varargin)
        if ( iscell(C{a}) )
            C{a} = C{a}{1};
        end 
        if ( isempty( C{a} ) )
            format = getFormat(varargin);
            error( 'Missing %s in %s specification. Correct format is: %s.', varargin{a}, upper(stage), format );
        end       
        if ( isnan( C{a} ) )
            format = getFormat(varargin);
            error( 'Missing numeric %s parameter in %s specification. Correct format is: %s.', varargin{a}, upper(stage), format );
        end        
    end
end

function str = getFormat(varargin)
    str = sprintf( '[%s], ', varargin{1}{:} );
    str = str(1:end-2);
end