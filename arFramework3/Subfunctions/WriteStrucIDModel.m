function [] = WriteStrucIDModel(modelname)
global ar
m = 1;
matVer = ver('MATLAB');
matlab_version = str2double(matVer.Version);

if ~exist('ar','var') && ~isfield(ar,'model')
    error('Please load and compile the model');
end

model_name = ar.model(m).name;

fprintf('\n Creating the %s model ... ', model_name);

p = ar.model(m).p(:);      % all parameters
fp = ar.model(m).fp(:);    % model conditions
x = ar.model(m).x;         % state variables
x0 = ar.model(m).px0;      % initial condition parameters
fx = ar.model(m).fx;       % initial condition equations
fy = ar.model(m).fy;       % observable equations
y = ar.model(m).y;         % observable names
z = ar.model(m).z;         % derived variables
fz = ar.model(m).fz;       % derived equations
u = ar.model(m).u;         % input variables
fu = ar.model(m).fu;       % input equation

fitted_pars = ar.pLabel(ar.qFit==1);

z = str2sym(z);
fz = str2sym(fz);

idx = find(~strcmp(p,fp));
p_cond = str2sym(p(idx));     % parameters in CONDITIONS
fp_cond = str2sym(fp(idx));


% output
h = str2sym(fy);
h = arSubs(h, z, fz, matlab_version);                 % substitute model conditions
h = arSubs(h, p_cond, fp_cond, matlab_version);       % substitute model conditions

% reaction equations
f = str2sym(fx);
f = arSubs(f, z, fz, matlab_version);                 % substitute drived equations
f = arSubs(f, p_cond, fp_cond, matlab_version);       % substitute model conditions


% known input
u = str2sym(u);
u = arSubs(u, z, fz, matlab_version);                 % substitute drived equations
u = arSubs(u, p_cond, fp_cond, matlab_version);       % substitute model conditions


% initial conditions
[q,idx] = ismember(x0,p);
if ~all(q==1)
    error('Error in initial conditions')
end
fx0 = fp(idx);

idx = ~ismember(x0,fx0);                             % replace initial conditions by their values in CONDITIONS

x0(idx) = fx0(idx);
x0 = str2sym(x0);
x0 = arSubs(x0, z, fz, matlab_version);              % substitute drived equations
x0 = arSubs(x0, p_cond, fp_cond, matlab_version);    % substitute model conditions


% parameters
pars = unique([symvar(h) symvar(f) symvar(x0) symvar(str2sym(setdiff(p, arFindPar(ar, {'sd'}, 'names'))))]);
pars = setdiff(pars, symvar(str2sym(x)));

h = arrayfun(@char, h, 'uniform', 0);
f = arrayfun(@char, f, 'uniform', 0);
u = arrayfun(@char, u, 'uniform', 0);
x0 = arrayfun(@char, x0, 'uniform', 0);
pars = arrayfun(@char, pars, 'uniform', 0);
%psyms = arrayfun(@char, psyms, 'uniform', 0);

% find fixed parameters
fix_par = ar.pLabel(ar.qFit==0);
fix_par_val = ar.p(ar.qFit==0);

fprintf('[done]\n');
%%
dfile = append(modelname, '.txt');
if exist(dfile, 'file') ; delete(dfile); end
diary(dfile)
diary on

fprintf('%s\n\n', 'Algebraic Rules!')
for entry = 1:length(fix_par)
    fprintf('%s = ', str2sym(fix_par(entry)))
    fprintf('%f\n', fix_par_val(entry))
end

fprintf('\n%s\n', 'ODEs (define the individual ODE equations - 1 per line)!')
for entry = 1:length(f)
    fprintf('d%s/dt = ', str2sym(x(entry)))
    fprintf('%s\n', str2sym(f(entry)))
end

fprintf('\n%s\n', 'Input variables!')
for entry = 1:length(u)
    fprintf('%s = ', str2sym(u(entry)))
    fprintf('%s\n', str2sym(fu(entry)))
end

fprintf('\n%s\n', 'Measured Outputs (define the measured sensors - 1 per line)!')
for entry = 1:length(y)
    fprintf('%s = ', str2sym(y(entry)))
    fprintf('%s\n', str2sym(fy(entry)))
end

fprintf('\n%s\n', 'Parameter names and values (define all the system parameters - 1 per line, OPTIONAL - define known paramter values)!')
for entry = 1:length(pars)
    fprintf('%s = ', str2sym(pars(entry)))
    fprintf('%s\n', '')
end
% for entry = 1:length(u)
%     fprintf('%s = \n', str2sym(fu(entry)))
% end

fprintf('\n%s\n', 'State names and initial values (define all the model state names - 1 per line, OPTIONAL - define known initial values)!')
for entry = 1:length(x)
    fprintf('%s = ', str2sym(x(entry)))
    fprintf('%s\n', str2sym(x0(entry)))
end

fprintf('\n%s\n', 'Analyse (list the unknown parameter and initial conditions which should be included into the strucutral identifiability analysis)!')
for entry = 1:length(fitted_pars)
    fprintf('%s\n', str2sym(fitted_pars(entry)))
end

diary off

end
