% out = arSubsRepeated(in, old, new, matlab_version)
% Better Symbolic substitution. Substitute until no more changes.
% (for self-substitutions of derived variables)
% 
% in    symbolic expression
% old   sym/str variable to be replaced
% new   sym/num variable with which to replace
%
% Example:
%       subs(cos(a)+sin(b),{a,b},{sym('alpha'),2}) returns
%       cos(alpha)+sin(2)
%
% See also ARSUBS, SUBS, SYM

function out = arSubsRepeated(in, old, new, matlab_version)

if(nargin<4)
    matVer = ver('MATLAB');
    matlab_version = str2double(matVer.Version);   
end

done = false;

old = arMyStr2Sym(old);
new = arMyStr2Sym(new);
in  = arMyStr2Sym(in);

k = 0; orig = in;
while ( ~done )
    out = arSubs(in, old, new, matlab_version);
    
    if ( k > 15 )
        v = '';
        for c = 1 : length( orig )
            if ~isequal( in(c), out(c) )
                v = sprintf( '%s\n%s = %s', v, char(old(c)), char(orig(c)) );
            end
        end
        error( 'Substitution recursion limit (15) exceeded!\nSolutions that cannot be obtained by simple substitution are not supported.\nDo you have any cyclic substitutions?\n%s\n', v );
    end
    
    % No more changes?
    if ( isempty( setdiff(out,in) ) )
        done = true;
    else
        in = out;
    end
    k = k + 1;
end

q = out;
out = cell(1,length(q));
for a = 1 : length(q)
    out{a} = char(q(a));
end
