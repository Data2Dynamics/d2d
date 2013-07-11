% Export to Mathematica
%
% arExportMathematica(m, d)
%
% ----------------------------------
% System: list of two lists.
% First, a list of derivatives of variables like {x1'[t] == k2*x2[t]+u1[t],
% x2'[t]==...}
% Second, a list of measured signals like {x2[t],x1[t] + k1*x3[t]^2/2, ...}
% 
% Variables: list of variable names like {x1,x2, ...}
% 
% Parameters: list of parameters like {k1,k2, ...}
% 
% U: If the system has inputs, then write them as a list of inputs like
% {u1[t],u2[t],
% ...}
% ----------------------------------
% johkar@fcc.chalmers.se


function arExportMathematica(m, d)

global ar

if(~exist([cd '/Mathematica' ], 'dir'))
    mkdir([cd '/Mathematica' ])
end


fid = fopen(sprintf('./Mathematica/%s.txt', ar.model(m).data(d).name), 'w');

c = ar.model(m).data(d).cLink;

fprintf(fid, 'states = {\n');
for jx=1:length(ar.model(m).x)
    strtemp = char(ar.model(m).condition(c).sym.fx(jx));
    for jv=1:length(ar.model(m).fv)
        strtemp = strrep(strtemp, sprintf('v[%i]', jv), char(ar.model(m).condition(c).sym.fv(jv)));
    end
    
    for jx2=1:length(ar.model(m).x)
        strtemp = strrep(strtemp, sprintf('x[%i]', jx2), sprintf('%s[%s]', strrep(ar.model(m).x{jx2}, '_', ''), ar.model(m).t));
    end
    
    for ju=1:length(ar.model(m).u)
        strtemp = strrep(strtemp, sprintf('u[%i]', ju), sprintf('%s[%s]', strrep(ar.model(m).u{ju}, '_', ''), ar.model(m).t));
    end
    
    for jp=1:length(ar.model(m).condition(c).p)
        strtemp = strrep(strtemp, sprintf('p[%i]', jp), sprintf('%s', strrep(ar.model(m).condition(c).p{jp}, '_', '')));
    end
    
    fprintf(fid, '\t%s''[%s] == %s,\n', strrep(ar.model(m).x{jx}, '_', ''), ar.model(m).t, strtemp);
end
fprintf(fid, '\t\n');

for jx=1:length(ar.model(m).x)    
    strtemp = char(ar.model(m).condition(c).sym.fpx0(jx));
    for jp=1:length(ar.model(m).condition(c).p)
        strtemp = strrep(strtemp, sprintf('p[%i]', jp), sprintf('%s', strrep(ar.model(m).condition(c).p{jp}, '_', '')));
    end
    for ju=1:length(ar.model(m).u)
        strtemp = strrep(strtemp, sprintf('u[%i]', ju), sprintf('%s[%s]', strrep(ar.model(m).u{ju}, '_', ''), ar.model(m).t));
    end
    fprintf(fid, '\t%s[0] == %s', strrep(ar.model(m).x{jx}, '_', ''), strtemp);
    
    if(jx<length(ar.model(m).x))
        fprintf(fid, ',\n');
    end
end
fprintf(fid, '\n}\n\n');

fprintf(fid, 'measurements = {\n');
for jy=1:length(ar.model(m).data(d).sym.fy)
    strtemp = char(ar.model(m).data(d).sym.fy(jy));
    
    for jx2=1:length(ar.model(m).x)
        strtemp = strrep(strtemp, sprintf('x[%i]', jx2), sprintf('%s[%s]', strrep(ar.model(m).x{jx2}, '_', ''), ar.model(m).t));
    end
    
    for ju=1:length(ar.model(m).u)
        strtemp = strrep(strtemp, sprintf('u[%i]', ju), sprintf('%s[%s]', strrep(ar.model(m).u{ju}, '_', ''), ar.model(m).t));
    end
    
    for jp=1:length(ar.model(m).data(d).p)
        strtemp = strrep(strtemp, sprintf('p[%i]', jp), sprintf('%s', strrep(ar.model(m).data(d).p{jp}, '_', '')));
    end
    
    fprintf(fid, '\t%s', strtemp);
    if(jy<length(ar.model(m).data(d).sym.fy))
        fprintf(fid, ',\n');
    end
end
fprintf(fid, '\n}\n\n');

fprintf(fid, 'variables = {\n');
for jx=1:length(ar.model(m).x)
    fprintf(fid, '\t%s', strrep(ar.model(m).x{jx}, '_', ''));
    if(jx<length(ar.model(m).x))
        fprintf(fid, ',\n');
    end
end
fprintf(fid, '\n}\n\n');

fprintf(fid, 'parameters = {\n');
ip = find(~ismember(ar.model(m).data(d).p, ar.model(m).data(d).pystd)); %R2013a compatible
for jp=ip
    fprintf(fid, '\t%s', strrep(ar.model(m).data(d).p{jp}, '_', ''));
    if(jp ~= ip(end))
        fprintf(fid, ',\n');
    end
end
fprintf(fid, '\n}\n\n');

fprintf(fid, 'inputs = {\n');
for ju=1:length(ar.model(m).u)
    fprintf(fid, '\t%s', strrep(ar.model(m).u{ju}, '_', ''));
    if(ju<length(ar.model(m).u))
        fprintf(fid, ',\n');
    end
end
fprintf(fid, '\n}\n\n');


fclose(fid);
