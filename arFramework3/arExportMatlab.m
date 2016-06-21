function arExportMatlab(m)

global ar

if(~exist('m','var'))
    m = 1;
end

%% write ode15s rhs file

fname = sprintf('%s_ode', ar.model(m).name);

fid = fopen([fname '.m'], 'w');

x = ar.model.x;
x0 = ar.model.px0;
u = ar.model.u;
z = ar.model.z;
p = setdiff(setdiff(ar.model.px, ar.model.px0), z);
f = ar.model.fx;

fpx0 = align_position(ar.model.fp, ar.model.p, x0);
fp = align_position(ar.model.fp, ar.model.p, p);

% subs u and z functions in rhs
fu = ar.model.fu;
fz = ar.model.fz;
f = subs(f, u, fu);
f = subs(f, z, fz);

fprintf(fid, 'function dxdt = %s(t,x,p,x0)\n\n', fname);
fprintf(fid, 'dxdt = zeros(size(x));\n\n');

% x
for j=1:length(x)
    fprintf(fid, '%s = x(%i);\n', x{j}, j);
end
fprintf(fid, '\n');

% p
for j=1:length(p)
    fprintf(fid, '%s = p(%i);\n', p{j}, j);
end
fprintf(fid, '\n');

% x0
for j=1:length(x0)
    fprintf(fid, '%s = x0(%i);\n', x0{j}, j);
end
fprintf(fid, '\n');

% f
for j=1:length(f)
    fprintf(fid, 'dxdt(%i) = %s;\n', j, char(f(j)));
end
fprintf(fid, '\n');

fclose(fid);


%% write ode15s driver script

fname2 = sprintf('%s_driver', ar.model(m).name);

fid = fopen([fname2 '.m'], 'w');

fprintf(fid, 'clc; clear;\n\n');

% p num
for j=find(ar.qDynamic)
    if(ar.qLog10(j))
        fprintf(fid, '%s = %g;\n', ar.pLabel{j}, 10^ar.p(j));
    else
        fprintf(fid, '%s = %g;\n', ar.pLabel{j}, ar.p(j));
    end
end
fprintf(fid, '\n');

% x0
for j=find(~cellfun(@strcmp, x0, fpx0))
    fprintf(fid, '%s = %s;\n', x0{j}, fpx0{j});
end
fprintf(fid, '\n');
fprintf(fid, 'x0 = [\n');
for j=1:length(x0)
    fprintf(fid, '\t%s\n', x0{j});
end
fprintf(fid, '\t];\n\n');

% p
for j=find(~cellfun(@strcmp, p', fp))
    fprintf(fid, '%s = %s;\n', p{j}, fp{j});
end
fprintf(fid, '\n');
fprintf(fid, 'p = [\n');
for j=1:length(p)
    fprintf(fid, '\t%s\n', p{j});
end
fprintf(fid, '\t];\n\n');

fprintf(fid, '[t,x] = ode15s(@(t,x) %s(t,x,p,x0), [%f %f], x0);\n\n', fname, ar.model.tLim(1), ar.model.tLim(2));

% state labels
fprintf(fid, 'state_label = {\n');
for j=1:length(ar.model(m).x)
    fprintf(fid, '\t%s\n', ar.model(m).x{j});
end
fprintf(fid, '\t};\n');

fprintf(fid, 'plot(t,x);\n\n');
fprintf(fid, 'legend(strrep(state_label, ''_'', ''\\_'');\n\n');

fclose(fid);


