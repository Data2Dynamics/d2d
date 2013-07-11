% Apply the quasi-steady-state assumtion after 
% Schauer and Heinrich (1982)
% Mathematical Biosciences
%
% function arQuasiSSApproximation(m)
%
% m:    model index

function arQuasiSSApproximation(m, fastreversibleindices, dependentvars)

global ar

qfastreversible = 1:length(ar.model(m).fv);
qfastreversible = ismember(qfastreversible, fastreversibleindices); %R2013a compatible

qdependentvars = 1:length(ar.model(m).x);
qdependentvars = ismember(qdependentvars, dependentvars); %R2013a compatible

independentvars = find(~qdependentvars);

n = length(ar.model(m).x);
r = length(dependentvars);

v = ar.model(m).fv(~qfastreversible);
w = ar.model(m).fv(qfastreversible);

R = ar.model(m).N(:,~qfastreversible);
S = ar.model(m).N(:,qfastreversible);

xs = mysym(ar.model(m).x);
ws = mysym(w);
dwdx = jacobian(ws,xs);

% check condition (26)
q = rank(S);
if(rank(dwdx)~=q)
    error('condition (26) is not fullfilled');
end 

rankSDwDx = double(rank(S*dwdx));
if(rankSDwDx ~= r)
    error('required number of dependent variables is %i', rankSDwDx);
end 

Ss = mysym(S);
Problem = Ss*ws;

ProblemStr = {};
for j=1:numel(Problem)
    if(Problem(j)~=0)
        ProblemStr{end+1} = char(Problem(j)); %#ok<AGROW>
    end
end
Solution = solve(ProblemStr{:}, ar.model(m).x{dependentvars});

% setup equilibrium relations G 
G = cell(r,1);
if(r>1)
    for j=1:r
        G{j} = getfield(Solution, ar.model(m).x{dependentvars(j)}); %#ok<GFLD>
    end
else
    G{1} = Solution;
end
Gs = mysym(G);
% replace equilibrium constants in G
for j=fastreversibleindices
    pname = ar.model(m).fv_ma_reverse_pbasename{j};
    Gs = mysubs(Gs, mysym([pname 'b']), mysym([pname 'f/' pname 'kD']));
end
fprintf('\nequilibrium relations:\n');
for j=1:r
    fprintf('\t%s\t= %s\n', ar.model(m).x{dependentvars(j)}, char(Gs(j)));
end

S1 = S(1:q,1:q);
% S2 = S(1:q,q+1:end);
S3 = S(q+1:end,1:q);
% S4 = S(q+1:end,q+1:end);

% linear independence of pool variables
if(rank(S1)~=size(S1,2))
    error('pool variables not linearly independent');
end 

% user choice
A2 = eye(n-q);
% A2 = [1 1 0; 0 1 0; 1 0 1];

A1 = - A2*S3*inv(S1); %#ok<*MINV>
A = [A1 A2];

Y = A * transpose(xs);
fprintf('\npooled variabled:\n');
for j=1:length(Y)
    fprintf('\t%s\n', char(Y(j)));
end

% differential equations for pooled variables
dgdxu = jacobian(Gs, xs(~qdependentvars));
T = [dgdxu;eye(n-r)]; % see (in parts) equation (37)

% replace equilibrium relations in v
vs = mysym(v);
vs = mysubs(vs, xs(dependentvars), G);
% replace equilibrium constants in vs
for j=fastreversibleindices
    pname = ar.model(m).fv_ma_reverse_pbasename{j};
    vs = mysubs(vs, mysym([pname 'b']), mysym([pname 'f/' pname 'kD']));
end

dxdtnew = inv(A*T) * A * R * vs;
% fprintf('simplifing expressions...\n');
% dxdtnew = simple(dxdtnew);
fprintf('\nnew ODE system:\n');
for j=1:length(dxdtnew)
    fprintf('\td/dt %s = %s\n', ar.model(m).x{independentvars(j)}, char(dxdtnew(j)));
end
fprintf('\nquantites:\n');
disp(findsym(dxdtnew));


% better sym
function out = mysym(in)
global ar
if(ar.config.isMaple) % maple
    out = sym(zeros(size(in)));
    for j1=1:size(in,1)
        for j2=1:size(in,2)
            out(j1,j2) = sym(in{j1,j2});
        end
    end
else % mupad
    out = sym(in);
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