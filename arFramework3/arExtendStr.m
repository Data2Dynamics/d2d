% Extends a string to a specified length
% 
% function str = arExtendStr(str, n, formatting)
%
%   str         - String to extend
%   n           - Number of characters to extend it to
%   formatting  - Text alignment (optional):
%                   'left' (default), 'right', 'center'
%
function str = arExtendStr(str, n, formatting)
    if ~exist( 'formatting', 'var' )
        formatting = 'left';
    end
    switch( formatting )
        case 'left'
            nd = n - length(str);
            S = ' ';
            str = [str S(ones(1,nd))];
        case 'right'
            nd = n - length(str);
            S = ' ';
            str = [S(ones(1,nd)) str];
        case 'center'
            nd = n - length(str);
            nl = floor( nd / 2 );
            nr = nd - nl;
            S = ' ';            
            str = [S(ones(1,nl)) str S(ones(1,nr))];
    end