% arExportMatlab(m)
% 
%   m       model index, default: 1
% 
%   This function writes two matlab files
%       1) one file representing the ODE system which can be called by
%       ode15s
%       2) a file for integrating the odes via ode15s
% 
% Examples:
%   arExportMatlab  % ar.model(1) is exported
% 
%   arExportMatlab(1:length(ar.model)) % all models are exported
% 

function arExportMatlab(m)

global ar

if(~exist('m','var'))
    m = 1;
end

if(length(m)>1)  % several models
    for i=1:length(m)
        arExportMatlab(m(i));
    end
    
else % a single model
    %% write ode15s rhs file
    
    fname = sprintf('%s_ode', ar.model(m).name);
    
    fid = fopen([fname '.m'], 'w');
    
    x = ar.model(m).x;
    x0 = ar.model(m).px0;
    u = ar.model(m).u;
    z = ar.model(m).z;
    p = setdiff(setdiff(union(ar.model(m).pu,ar.model(m).px), ar.model(m).px0), z);
    f = ar.model(m).fx;
    
    fpx0 = align_position(ar.model(m).fp, ar.model(m).p, x0);
    fp = align_position(ar.model(m).fp, ar.model(m).p, p);
    
    % subs u and z functions in rhs
    fu = ar.model(m).fu;
    fz = ar.model(m).fz;
    if(size(z,1)==1 && size(fz,1)>1)
        fz = fz';
    end
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
    
    fprintf(fid, '[t,x] = ode15s(@(t,x) %s(t,x,p,x0), [%f %f], x0);\n\n', fname, ar.model(m).tLim(1), ar.model(m).tLim(2));
    
    % state labels
    fprintf(fid, 'state_label = {\n');
    for j=1:length(ar.model(m).x)
        fprintf(fid, '\t%s\n', ar.model(m).x{j});
    end
    fprintf(fid, '\t};\n');
    
    fprintf(fid, 'plot(t,x);\n\n');
    fprintf(fid, 'legend(strrep(state_label, ''_'', ''\\_''));\n\n');
    
    fclose(fid);
        
end


function compile_arInputFunctionsC

