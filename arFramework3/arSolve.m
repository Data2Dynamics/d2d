% solve ODE system analytically
%
% arSolve(m, c)
%   m:              model index         
%   c:              condition index

function arSolve(m, c)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

rhs = (ar.model(m).N.*ar.model(m).condition(c).sym.C)*ar.model(m).condition(c).sym.fv;
rhs = arSubs(rhs, ar.model(m).us, ar.model(m).condition(c).sym.fu);

fstrings = cell(1, 2*length(ar.model(m).x));
for j=1:length(ar.model(m).x)
    fstrings{j} = sprintf('D%s = %s', ar.model(m).xs{j}, char(rhs(j)));
end
for j=1:length(ar.model(m).xs)
    fstrings{j+length(ar.model(m).xs)} = sprintf('%s(0) = %s', ar.model(m).xs{j}, ...
        char(ar.model(m).condition(c).sym.fpx0(j)));
end
fstrings = strrep(fstrings, '[', '');
fstrings = strrep(fstrings, ']', '');

for j=1:length(fstrings)
    disp(fstrings{j})
end
ar.model(m).condition(c).xt_solution = dsolve(fstrings{:});

if(~isempty(ar.model(m).condition(c).xt_solution))
    fprintf('---\n');
    for j=1:length(ar.model(m).x)
        xtmp = char(ar.model(m).sym.xs(j));
        xtmp = strrep(xtmp, '[', '');
        xtmp = strrep(xtmp, ']', '');
        tmpsym = getfield(ar.model(m).condition(c).xt_solution, xtmp); %#ok<GFLD>
        fprintf('%s(t) = %s\n', xtmp, char(tmpsym));
    end
end


