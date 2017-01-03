% solve steady state condition
%
% arSolveSS(m, fu, fextra)
%
% m:        model index
% fu:       custom initial input
% fextra:   additional equations

function arSolveSS(m, fu, fextra)

global ar

fprintf('solving steady state condition...\n');

rhs = sym(ar.model(m).fx);
if(exist('fu','var'))
    rhs = arSubs(rhs, ar.model(m).sym.u, sym(fu));
elseif(~isempty(ar.model(m).fu))
    rhs = arSubs(rhs, ar.model(m).sym.u, sym(ar.model(m).fu));
end
rhs = arSubs(rhs, ar.model(m).sym.x, ar.model(m).sym.px0);

qnonzero = rhs~=0;
fstrings = cell(1, sum(logical(qnonzero)));
dxdt = cell(1, sum(logical(qnonzero)));
fcount = 1;
for j=find(qnonzero)'
    fstrings{fcount} = sprintf('%s = 0', char(rhs(j)));
    dxdt{fcount} = ar.model(m).x{j};
    fcount = fcount + 1;
end

fp = cell(1, sum(logical(qnonzero)));
fcount = 1;
for j=find(qnonzero)'
    fp{fcount} = ar.model(m).px0{j};
    fcount = fcount + 1;
end

fprintf('\nsteady state equations (%i):\n', length(fstrings));
for j=1:length(fstrings)
    disp(['d' dxdt{j} '/dt = ' fstrings{j}])
end

fprintf('\nsolving for (%i):\n', length(fp));
for j=1:length(fp)
    disp(fp{j})
end
fprintf('\n---\n');

if(exist('fextra','var'))
    ar.model(m).ss_solution = solve(fstrings{:}, fextra{:}, fp{:});
else
    ar.model(m).ss_solution = solve(fstrings{:}, fp{:});
end

if(~isempty(ar.model(m).ss_solution))
    fprintf('---\n');
    tmpstr = getfield(ar.model(m).ss_solution, fp{1}); %#ok<GFLD>
    
    for jj = 1:length(tmpstr)
        fprintf('\n%s Solution %i:\n', ar.config.comment_string, jj);
        for j = 1:length(fp)
            tmpstr = getfield(ar.model(m).ss_solution, fp{j}); %#ok<GFLD>
            %fprintf('%s\t  "%s"\n', fp{j}, char(simple(tmpstr(jj))));
            fprintf('%s\t  "%s"\n', fp{j}, char(tmpstr(jj)));
        end
    end
end

