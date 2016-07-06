% arExportMatlab(m)
% 
%   m       model index, default: 1
% 
%   This function writes two matlab files
%       1) one file representing the ODE system which can be called by
%       ode15s
%       2) a file for integrating the odes via ode15s
% 
% The function has been tested for Raia, Becker, Bachmann, Boehm.
% Swameye does not work because the c-function for splines occuring in the
% ODEs is not available.
% 
% Examples:
%   arExportMatlab  % ar.model(1) is exported
% 
%   arExportMatlab(1:length(ar.model)) % all models are exported
% 

function arExportMatlab(m)


global ar

if(~exist('m','var') || isempty(m))
    m = 1;
end

if(length(m)>1)  % several models
    for i=1:length(m)
        arExportMatlab(m(i));
    end
    
else % a single model
    if strcmp(version('-release'), '2016a')
        warning('off', 'symbolic:sym:sym:DeprecateExpressions')
    end
    %% write ode15s rhs file
    
    fname = sprintf('%s_ode', ar.model(m).name);
    
    fid = fopen([fname '.m'], 'w');
    
    x = ar.model(m).x;
    x0 = ar.model(m).px0;
    u = ar.model(m).u;
    z = ar.model(m).z;
    p = setdiff(setdiff(union(ar.model(m).pu,ar.model(m).px), ar.model(m).px0), z);
    f = ar.model(m).fx;
        
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
    
    fprintf(fid, 'clc; clear; close all; \n\n');
    
    % p num
    fprintf(fid,'%s parameters and/or initial conditions:\n','%%');
    fprintf(fid,'p = struct;\n');
    for j=find(ar.qDynamic)
        if(ar.qLog10(j))
            fprintf(fid, 'p.%s = %g;\n', ar.pLabel{j}, 10^ar.p(j));  % putting the model parameters in a struct is required to prevent overwriting and multiple replacements within conditions
        else
            fprintf(fid, 'p.%s = %g;\n', ar.pLabel{j}, ar.p(j));
        end
    end
    fprintf(fid, '\n');
    
    fprintf(fid, '%s loop over all "conditions", i.e. over all settings with different initial values and/or over pertubations \n','%%');
    fprintf(fid, 'for jc=1:%i\n\t \n',length(ar.model(m).condition));
    fprintf(fid,'\tswitch jc\n');
    
    for c = 1:length(ar.model(m).condition)
        xda = ~isempty(ar.model(m).x);
        uda = ~isempty(ar.model(m).u);
        
        c_p = ar.model(m).condition(c).pold;
        c_fp = ar.model(m).condition(c).fp;
        c_fp = ReplaceVarnames(c_fp);  % substitute each variable, e.g. x by p.x
        
        u = ar.model(m).condition(c).fu;
        
        fprintf(fid, '\t\tcase %i  %s condition/pertubation %i\n',c,'%',c);
        fprintf(fid,'\t\t\t%s condition-specific parameters:\n','%');

        for i=1:length(c_p)
%             if(strcmp(c_p{i},c_fp{i})~=1)
                fprintf(fid,'\t\t\t%s = %s;\n',c_p{i},c_fp{i});
%             end
        end

        if xda
            fprintf(fid, '\n');
            fprintf(fid, '\t\tx0 = [');
            for j=1:length(x)
                fprintf(fid, 'init_%s ', x{j});
                if rem(j,10)==0
                    fprintf(fid,'...\n\t\t\t');
                else
                    fprintf(fid,' ');
                end
            end
            fprintf(fid, '];\n\n');
        end
%         
%         % p
%         for j=find(~cellfun(@strcmp, p', fp))
%             fprintf(fid, '%s = %s;\n', p{j}, fp{j});
%         end
%         fprintf(fid, '\n');
        fprintf(fid, '\t\tpcond = [');
        for j=1:length(p)
            fprintf(fid, '%s', p{j});
            if rem(j,10)==0
                fprintf(fid,'...\n\t\t\t');
            else
                fprintf(fid,' ');
            end
        end
        fprintf(fid, '];\n\n');
        
        ds = ar.model(m).condition(c).dLink;
        pred_name = cell(0);
%         preds = [];
        for d=ds
            pred_name{end+1} = ar.model(m).data(d).t;
%             preds = union(preds,ar.model(m).data(d).tFine);
        end
        pred_name = unique(pred_name);
        if length(pred_name)>1
            pred_name
            error('arExportMatlab currently handles only one predictor per model')
        else       
            pred_name = pred_name{1};
        end
                
        fprintf(fid, '\t\t%s = [',pred_name);
%         preds = unique(preds);
        preds = ar.model(m).condition(c).tFine; % attention: prediciton values (horizontal axis) can be condition-dependent
        for it=1:length(preds)
            fprintf(fid,'%f',preds(it));
            if rem(it,10)==0
                fprintf(fid,'...\n\t\t\t');
            else
                fprintf(fid,' ');
            end
        end
        fprintf(fid,'];\n');

        
        if xda
%             fprintf(fid, '\t\t[t,x] = ode15s(@(t,x) %s(t,x,p,x0), [%f %f], x0);\n\n', fname, ar.model(m).tLim(1), ar.model(m).tLim(2));
            fprintf(fid, '\n\t\t[t,x] = ode15s(@(t,x) %s(t,x,pcond,x0), %s, x0);\n', fname, pred_name);
        else
            fprintf(fid, '\t\tx=[];\n');
        end
        if uda            
            fprintf(fid, '\n');
            fprintf(fid, '\t\tu = [');
            for j=1:length(u)
                fprintf(fid,'%s  ',strrep(strrep(strrep(u{j},'^','.^'),'/','./'),'*','.*'));
            end
            fprintf(fid, '];\n');

        else
            fprintf(fid, '\t\tu=[];\n');            
        end
    end
    
    fprintf(fid,'\n\tend %s switch over conditions\n','%');
    
    
    % state labels
    fprintf(fid, '\tstate_labels = {\n');
    for j=1:length(ar.model(m).x)
        fprintf(fid, '\t''%s''\n', ar.model(m).x{j});
    end
    fprintf(fid, '\t\t};\n');
    
    fprintf(fid, '%s plotting \n','%%');
    fprintf(fid, '\tif(~isempty(x))\n');
    fprintf(fid, '\t\tfigure;\n');
    fprintf(fid, '\t\tplot(t,x);\n');
    fprintf(fid, '\t\tlegend(strrep(state_labels, ''_'', ''\\_''));\n');
    fprintf(fid, '\tend\n\n');

    fprintf(fid, '\tu_labels = {\n');
    for j=1:length(ar.model(m).u)
        fprintf(fid, '\t''%s''\n', ar.model(m).u{j});
    end
    fprintf(fid, '\t\t};\n');

    fprintf(fid, '\tif(~isempty(u))\n');    
    fprintf(fid, '\tif length(u)==1, u=u*ones(size(%s));end\n',pred_name);
    fprintf(fid, '\t\tfigure;\n');
    fprintf(fid, '\t\tplot(%s,u);\n',pred_name);
    fprintf(fid, '\t\tlegend(strrep(u_labels, ''_'', ''\\_''));\n');
    fprintf(fid, '\t\txlabel(''%s'');\n',strrep(pred_name,'_','\_'));    
    fprintf(fid, '\tend\n\n');
    
    fprintf(fid,'end\n');
    fclose(fid);
        
    
    if strcmp(version('-release'), '2016a')
        warning('on', 'symbolic:sym:sym:DeprecateExpressions')
    end
    
end

function fp = ReplaceVarnames(fp)
if(iscell(fp))
    for i=1:length(fp)
        fp{i} = ReplaceVarnames(fp{i});
    end
else
    fp = sym(fp);
    vs = symvar(fp);
    if(~isempty(vs))
        for i=1:length(vs)
            fp = subs(fp,vs(i),['pppp_',char(vs(i))]);
        end
        fp = strrep(char(fp),'pppp_','p.');
    else
        fp = char(fp);
    end
end

function compile_arInputFunctionsC

