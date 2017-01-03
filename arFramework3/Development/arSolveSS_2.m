% solve steady state condition
%
% arSolveSS_2(m, fu)
%
% m:    model index
% fu:   custom initial input

function arSolveSS_2(m, fu)

global ar

fprintf('solving steady state condition...\n');

rhs = sym(ar.model(m).fx);
if(exist('fu','var'))
    rhs = arSubs(rhs, ar.model(m).sym.u, sym(fu));
elseif(~isempty(ar.model(m).fu))
    rhs = arSubs(rhs, sym(ar.model(m).u), zeros(1,length(ar.model(m).u)));
end
rhs = arSubs(rhs, sym(ar.model(m).x), sym(ar.model(m).px0));

rhs= arSubs(rhs, sym(ar.pLabel(strncmp('init_',ar.pLabel,5) & ar.qFit==2 & ar.p == 0)), zeros(1,length(sym(ar.pLabel(strncmp('init_',ar.pLabel,5) & ar.qFit==2 & ar.p == 0)))));

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

tmp_pars = ar.pLabel';
tmp_list = cell(length(ar.pLabel),2);
tmp_list(:,2) = tmp_pars;
tmp_list(:,1) = num2cell((1:length(ar.pLabel))');
tmp_list
prompt = ['Press Enter to solve eq. for initial values or supply variable numbers as a vector of length %i from this list, e.g. [1 5 7] \n',sum(logical(qnonzero))];
input_pars = input(prompt);
if(isempty(input_pars))
    input_pars = find(strncmp(ar.pLabel,'init_',5) & ar.qFit~=2);
end
for j=1:length(input_pars)
    fp{fcount} = ar.pLabel{input_pars(j)};
    fcount = fcount + 1;
end

fprintf('\nsteady state equations (%i):\n', length(fstrings));
for j=1:length(fstrings)
    disp(fstrings{j})
end

fp = fp(~cellfun('isempty',fp));
fprintf('\nsolving for (%i):\n', length(fp));
for j=1:length(fp)
    disp(fp{j})
end
fprintf('\n---\n');

ar.model(m).ss_solution = solve(fstrings{:}, fp{:});

if(~isempty(ar.model(m).ss_solution))
    fprintf('---\n');
    tmpstr = getfield(ar.model(m).ss_solution, fp{1}); %#ok<GFLD>
    
    for jj = 1:length(tmpstr)
        fprintf('\n%s Solution %i:\n', ar.config.comment_string, jj);
        for j = 1:length(fp)
            tmpstr = getfield(ar.model(m).ss_solution, fp{j}); %#ok<GFLD>
            fprintf('%s\t  "%s"\n', fp{j}, char(simplify(tmpstr(jj))));
        end
    end
end


