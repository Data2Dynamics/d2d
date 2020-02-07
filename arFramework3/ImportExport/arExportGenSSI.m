% arExportGenSSI(m, d)
%
% Export to GenSSI identifiability software
%
%   m - Model index
%   d - Data index (you can find this index with arFindData) (ar.model(m).data(d))
% 
% arExportGenSSI(1,1)

function arExportGenSSI(m, d)

global ar

if(~exist([cd '/GenSSI' ], 'dir'))
    mkdir([cd '/GenSSI' ])
end


fid = fopen(sprintf('./GenSSI/%s.m', ar.model(m).data(d).name), 'w');

c = ar.model(m).data(d).cLink;

fprintf(fid, 'clc\n');
fprintf(fid, 'clear\n\n');

fprintf(fid, 'example_folder=''%s'';\n', ar.model(m).data(d).name);
fprintf(fid, 'result_folder=''%s'';\n\n', ar.model(m).data(d).name);

fprintf(fid, 'syms ');
for jx=1:length(ar.model(m).x)
    fprintf(fid, '%s ',  strrep(ar.model(m).x{jx}, '_', ''));
end
ip = find(~ismember(ar.model(m).data(d).p, ar.model(m).data(d).pystd)); %R2013a compatible
for jp=ip
    fprintf(fid, '%s ', strrep(ar.model(m).data(d).p{jp}, '_', ''));
end
fprintf(fid, '\n\n');

%% Number of derivatives
fprintf(fid, 'Nder=5;\n\n');

%% Number of states
fprintf(fid, 'Neq=%i;\n\n', length(ar.model(m).x));

%% Number of controls
fprintf(fid, 'Noc=0;\n\n');

%% Number of observables
fprintf(fid, 'Nobs=%i;\n\n', length(ar.model(m).data(d).y));

%% Number of model parameters
fprintf(fid, 'Npar=%i;\n\n', length(ip));

%% States
fprintf(fid, 'X = [');
for jx=1:length(ar.model(m).x)
    fprintf(fid, '%s ',  strrep(ar.model(m).x{jx}, '_', ''));
end
fprintf(fid, '];\n\n');

%% Equations of the model
for jx=1:length(ar.model(m).x)
    strtemp = char(ar.model(m).condition(c).sym.fx(jx));
    for jv=1:length(ar.model(m).fv)
        strtemp = strrep(strtemp, sprintf('v[%i]', jv), char(ar.model(m).condition(c).sym.fv(jv)));
    end
    
    for jx2=1:length(ar.model(m).x)
        strtemp = strrep(strtemp, sprintf('x[%i]', jx2), sprintf('%s', strrep(ar.model(m).x{jx2}, '_', '')));
    end
    
    for ju=1:length(ar.model(m).u)
        strtemp = strrep(strtemp, sprintf('u[%i]', ju), sprintf('%s', strrep(ar.model(m).u{ju}, '_', '')));
    end
    
    for jp=1:length(ar.model(m).condition(c).p)
        strtemp = strrep(strtemp, sprintf('p[%i]', jp), sprintf('%s', strrep(ar.model(m).condition(c).p{jp}, '_', '')));
    end
    
    fprintf(fid, 'A%i = %s;\n', jx, strtemp);
end
fprintf(fid, 'F = [');
for jx=1:length(ar.model(m).x)
    fprintf(fid, 'A%i ', jx);
end
fprintf(fid, '];\n\n');

%% Controls
fprintf(fid, 'G = [');
for jx=1:length(ar.model(m).x)
    fprintf(fid, '0 ', jx);
end
fprintf(fid, '];\n\n');

%% Observables
for jy=1:length(ar.model(m).data(d).sym.fy)
    strtemp = char(ar.model(m).data(d).sym.fy(jy));
    
    for jx2=1:length(ar.model(m).x)
        strtemp = strrep(strtemp, sprintf('x[%i]', jx2), sprintf('%s', strrep(ar.model(m).x{jx2}, '_', '')));
    end
    
    for ju=1:length(ar.model(m).u)
        strtemp = strrep(strtemp, sprintf('u[%i]', ju), sprintf('%s', strrep(ar.model(m).u{ju}, '_', '')));
    end
    
    for jp=1:length(ar.model(m).data(d).p)
        strtemp = strrep(strtemp, sprintf('p[%i]', jp), sprintf('%s', strrep(ar.model(m).data(d).p{jp}, '_', '')));
    end
    
    fprintf(fid, 'h%i = %s;\n', jy, strtemp);
end
fprintf(fid, 'H = [');
for jy=1:length(ar.model(m).data(d).sym.fy)
    fprintf(fid, 'h%i ', jy);
end
fprintf(fid, '];\n\n');

%% Initial conditions
fprintf(fid, 'IC = [');
for jx=1:length(ar.model(m).x)    
    strtemp = char(ar.model(m).condition(c).sym.fpx0(jx));
    for jp=1:length(ar.model(m).condition(c).p)
        strtemp = strrep(strtemp, sprintf('p[%i]', jp), sprintf('%s', strrep(ar.model(m).condition(c).p{jp}, '_', '')));
    end
    for ju=1:length(ar.model(m).u)
        strtemp = strrep(strtemp, sprintf('u[%i]', ju), sprintf('%s', strrep(ar.model(m).u{ju}, '_', '')));
    end
    fprintf(fid, '%s ', strtemp);
end
fprintf(fid, '];\n\n');

%% parameter for identifiability
fprintf(fid, 'Par = [');
for jp=ip
    fprintf(fid, '%s', strrep(ar.model(m).data(d).p{jp}, '_', ''));
    if(jp~=ip(end))
        fprintf(fid, ', ');
    end
end
fprintf(fid, '];\n\n');

%% call GenSSI
fprintf(fid, 'GenSSI_generating_series(F,G,Noc,Neq,Nder,Nobs,H,X,IC,Par,example_folder,result_folder);\n');

fclose(fid);
