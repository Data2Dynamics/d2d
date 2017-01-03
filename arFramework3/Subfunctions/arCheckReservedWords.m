% Checks whether a word reserved by the symbolic toolbox exists in a string

function arCheckReservedWords( strings, useLocation, offendingString )

    keywords = {'time','gamma','sin','cos','tan','beta','log','asin','atan','acos','acot','cot','theta','D','x'};
    inter = intersect(strings, keywords);
    if(~isempty(inter))
        if ( nargin > 2 )
            error('Fatal error: used reserved word(s) "%s" in %s (in %s).\nReserved words are: %s\n', sprintf('%s ', inter{:}), useLocation, offendingString, sprintf('%s ', keywords{:}) );
        else
            error('Fatal error: used reserved word(s) "%s" in %s.\nReserved words are: %s\n', sprintf('%s ', inter{:}), useLocation, sprintf('%s ', keywords{:}));
        end
    end

end
