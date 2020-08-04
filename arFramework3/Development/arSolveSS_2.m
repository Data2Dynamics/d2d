% solve steady state condition
%
% arSolveSS_2(m, [fu])
%
% m:    model index
% fu:   custom initial input

function arSolveSS_2(m, fu)
persistent ver % checking the version every time is extremely time-consuming.
if isempty(ver)
    ver = version;
    ver = str2num(ver(1:3));
end

global ar

fprintf('solving steady state condition...\n');

rhs = arSym(ar.model(m).fx);
if(exist('fu','var'))
    rhs = arSubs(rhs, ar.model(m).sym.u, arSym(fu));
elseif(~isempty(ar.model(m).fu))
    rhs = arSubs(rhs, arSym(ar.model(m).u), zeros(1,length(ar.model(m).u)));
end
rhs = arSubs(rhs, arSym(ar.model(m).x), arSym(ar.model(m).px0));

rhs= arSubs(rhs, arSym(ar.pLabel(strncmp('init_',ar.pLabel,5) & ar.qFit==2 & ar.p == 0)), zeros(1,length(arSym(ar.pLabel(strncmp('init_',ar.pLabel,5) & ar.qFit==2 & ar.p == 0)))));

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
prompt = sprintf(' Press Enter to solve equations for initial values OR\n supply variable numbers as a vector of length %i from this list, e.g. [1 5 7] \n',sum(logical(qnonzero)));
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

if ver<9.4  % new symbolic toolbox, I used the same distinction like in arSym without checking
    ar.model(m).ss_solution = solve(fstrings{:}, fp{:});
else
    fstrings = arSym(fstrings);
    fp = arSym(fp);
    ar.model(m).ss_solution = solve(fstrings, fp);
end

if(~isempty(ar.model(m).ss_solution))
    fprintf('---\n');
    for j = 1:length(fp)
        if ver<9.4  % new symbolic toolbox, I used the same distinction like in arSym without checking
            solu = getfield(ar.model(m).ss_solution, fp{j}); %#ok<GFLD>
        else
            solu = getfield(ar.model(m).ss_solution, char(fp(j))); %#ok<GFLD>
        end
        
        for jj = 1:length(solu)
            fprintf('\n%s Solution %i:\n', ar.config.comment_string, jj);
            if ver<9.4  % new symbolic toolbox, I used the same distinction like in arSym without checking
                fprintf('%s\t  "%s"\n', fp{j}, char(simplify(solu(jj))));
            else
                fprintf('%s\t  "%s"\n', char(fp(j)), char(simplify(solu(jj))));
            end
        end
    end
end


