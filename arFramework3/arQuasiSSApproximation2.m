% Apply the quasi-steady-state assumtion after 
% Stiefenhofer (1998)
% Journal of Mathmatical Biology
%
% function arQuasiSSApproximation2(m)
%
% m:    model index

function arQuasiSSApproximation2(m, fastreversibleindices, dependentvars)

global ar

qfastreversible = 1:length(ar.model(m).fv);
qfastreversible = ismember(qfastreversible, fastreversibleindices); %R2013a compatible

qdependentvars = 1:length(ar.model(m).x);
qdependentvars = ismember(qdependentvars, dependentvars); %R2013a compatible

independentvars = find(~qdependentvars);

R1 = ar.model(m).N(~qdependentvars, ~qfastreversible);
R2 = ar.model(m).N(qdependentvars, ~qfastreversible);

S1 = ar.model(m).N(~qdependentvars, qfastreversible);
S2 = ar.model(m).N(qdependentvars, qfastreversible);

v = sym(ar.model(m).fv(~qfastreversible));
w = sym(ar.model(m).fv(qfastreversible));

% calculate y(x) see (3.6)
ydotfast = S2*w;
ProblemStr = {};
for j=1:numel(ydotfast)
    if(ydotfast(j)~=0)
        ProblemStr{end+1} = char(ydotfast(j)); %#ok<AGROW>
    end
end
Solution = solve(ProblemStr{:}, ar.model(m).x{dependentvars});
yx = cell(length(dependentvars),1);
if(length(dependentvars)>1)
    for j=1:length(dependentvars)
        yx{j} = getfield(Solution, ar.model(m).x{dependentvars(j)}); %#ok<GFLD>
    end
else
    yx{1} = Solution;
end
fprintf('\ny(x):\n');
for j=1:length(dependentvars)
    fprintf('\t%s\t= %s\n', ar.model(m).x{dependentvars(j)}, char(yx{j}));
end

% calculate f(x,y(x),0) see fe[.] in (3.8)
xdotslow = R1*v;
for j=1:length(dependentvars)
    xdotslow = mysubs(xdotslow, ar.model(m).x{dependentvars(j)}, yx{j});
end
fprintf('\nxdotslow:\n');
for j=1:length(xdotslow)
    fprintf('\td[%s]/dt\t= %s\n', ar.model(m).x{independentvars(j)}, char(xdotslow(j)));
end



% better subs
function out = mysubs(in, old, new, flag)
if(~exist('flag','var'))
    flag = 0;
end
if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
    out = subs(in, old, new, flag);
else
    out = in;
end