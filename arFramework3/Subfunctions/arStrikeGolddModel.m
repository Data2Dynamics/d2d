
% this function generates the model which is neede by StrIkE-GOLDD(v3.0) toolbox for the structural identifiability
% and observability analysis.


function arStrikeGolddModel(m)

global ar

matVer = ver('MATLAB');
matlab_version = str2double(matVer.Version);

if ~exist('ar','var') && ~isfield(ar,'model')
    error('Please load and compile the model');
end

model_name = ar.model(m).name;
checksum = addToCheckSum(ar.model(m).name);

fprintf('\n Creating the %s model ... ', model_name);


% known and unknown inputs
if isfield(ar.model,'w') && ~isempty(ar.model(m).w)
    w = ar.model(m).w;                  % unknown inputs should be specified by the used in ar.model.w 
    idx = ~ismember(ar.model(m).u ,w);
    u = ar.model(m).u(idx);             % exlude unknown inputs from known inputs
else
    w = '';
    u = ar.model(m).u;
end


p = ar.model(m).p(:);      % all parameters
fp = ar.model(m).fp(:);    % model conditions
x = ar.model(m).x;         % state variables
x0 = ar.model(m).px0;      % initial condition parameters
fx = ar.model(m).fx;       % initial condition equations
fy = ar.model(m).fy;       % observable equations
z = ar.model(m).z;         % derived variables
fz = ar.model(m).fz;       % derived equations

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

% unknown input
w = str2sym(w);
w = arSubs(w, z, fz, matlab_version);                 % substitute drived equations
w = arSubs(w, p_cond, fp_cond, matlab_version);       % substitute model conditions



% initial conditions
[q,idx] = ismember(x0,p);
if ~all(q==1)
    error('Error in initial conditions (arStrikeGolddModel.m)')
end
fx0 = fp(idx);

idx = ~ismember(x0,fx0);                             % replace initial conditions by their values in CONDITIONS

x0(idx) = fx0(idx);
x0 = str2sym(x0);
x0 = arSubs(x0, z, fz, matlab_version);              % substitute drived equations
x0 = arSubs(x0, p_cond, fp_cond, matlab_version);    % substitute model conditions

% find known initial conditions
known_ics = cell(1,length(x0));
known_ics(:) = {'0'};
for i = 1:length(x0)
    if isempty(symvar(x0(i)))
            known_ics(i)={'1'};
    end
end


% parameters
pars = unique([symvar(h) symvar(f) symvar(x0)]);
pars = setdiff(pars, symvar(str2sym(x)));

% syms
psyms = unique([pars symvar(u) symvar(w) symvar(str2sym(p)) symvar(str2sym(x))]);


h = arrayfun(@char, h, 'uniform', 0);
f = arrayfun(@char, f, 'uniform', 0);
u = arrayfun(@char, u, 'uniform', 0);
w = arrayfun(@char, w, 'uniform', 0);
x0 = arrayfun(@char, x0, 'uniform', 0);
pars = arrayfun(@char, pars, 'uniform', 0);
psyms = arrayfun(@char, psyms, 'uniform', 0);



% add to check sum
checksum = addToCheckSum(psyms, checksum);
checksum = addToCheckSum(x, checksum);
checksum = addToCheckSum(h, checksum);
checksum = addToCheckSum(ar.model(m).u, checksum);
checksum = addToCheckSum(w, checksum);
checksum = addToCheckSum(pars, checksum);
checksum = addToCheckSum(f, checksum);
checksum = addToCheckSum(x0, checksum);
checksum = addToCheckSum(known_ics, checksum);


checkstr = getCheckStr(checksum);
ar.ia.checkstr = checkstr;

ar.ia.modelname = model_name;
ar.ia.modelname_checkstr = [model_name '_' checkstr];

ar.ia.data.syms = psyms;
ar.ia.data.x = x;
ar.ia.data.h = h;
ar.ia.data.u = u;
ar.ia.data.w = w;
ar.ia.data.p = pars;
ar.ia.data.f = f;
ar.ia.data.ics = x0;
ar.ia.data.known_ics = known_ics;

fprintf('[done]\n');

end


function checksum = addToCheckSum(str, checksum)
algs = {'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if(nargin<2)
    checksum = java.security.MessageDigest.getInstance(algs{2});
end

if iscell(str)
    str = [str{:}];
end
if(~isempty(str))
    checksum.update(uint8(str(:)));
end

end

function checkstr = getCheckStr(checksum)
h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

clear checksum
end