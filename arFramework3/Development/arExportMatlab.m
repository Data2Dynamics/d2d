% arExportMatlab([m])
% 
% This function writes two matlab files
%   1) one file representing the ODE system which can be called by ode15s
%   2) a file for integrating the odes via ode15s
% 
%   m   model index [1]
% 
% The function has been tested for Raia, Becker, Bachmann, Boehm.
% Swameye does not work because the c-function for splines occuring in the
% ODEs is not available.
% 
% Examples:
%   1) arExportMatlab  % ar.model(1) is exported
%   2) arExportMatlab(1:length(ar.model)) % all models are exported

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
    for j=1:length(ar.p)
        if(ar.qLog10(j))
            fprintf(fid, 'p_struct.%s = %g;\n', ar.pLabel{j}, 10^ar.p(j));  % putting the model parameters in a struct is required to prevent overwriting and multiple replacements within conditions
        else
            fprintf(fid, 'p_struct.%s = %g;\n', ar.pLabel{j}, ar.p(j));
        end
    end
    fprintf(fid, '\n');

    fprintf(fid,'p_labels = fieldnames(p_struct);\n');
    fprintf(fid,'p_cell = cell(size(p_labels));\n');

    %% predictors (normalerweise die Zeit)
    fprintf(fid,'\t%% predictors (often the time)\n');
    
    for d = 1:length(ar.model(m).data)
        pred_labels{d} = ar.model(m).data(d).t;
        fprintf(fid,'pred_struct.%s = [',pred_labels{d});                
        preds = ar.model(m).condition(ar.model(m).data(d).cLink).tFine; % attention: prediciton values (horizontal axis) can be condition-dependent
        for it=1:length(preds)
            fprintf(fid,'%f',preds(it));
            if rem(it,10)==0
                fprintf(fid,'...\n\t');
            else
                fprintf(fid,' ');
            end
        end
        fprintf(fid,'];\n');
    end
    fprintf(fid, '\npred_labels = fieldnames(pred_struct);\n');
    fprintf(fid, 'for i=1:length(pred_labels)\n');
    fprintf(fid, '\tpred_cell{i} = pred_struct.(pred_labels{i});\n');
    fprintf(fid, 'end\n');

    
    %% conditions
    fprintf(fid,'for i=1:length(p_labels)\n\tp_cell{i} = p_struct.(p_labels{i});\nend\n');
    
    fprintf(fid, '%s loop over all "conditions", i.e. over all settings with different initial values and/or over pertubations \n','%%');
    fprintf(fid, 'for jc=1:%i\n\t \n',length(ar.model(m).condition));    
    % state labels
    fprintf(fid, '\tx_labels = {\n');
    for j=1:length(ar.model(m).x)
        fprintf(fid, '\t''%s''\n', ar.model(m).x{j});
    end
    fprintf(fid, '\t\t};\n');

    fprintf(fid, '\tu_labels = {\n');
    for j=1:length(ar.model(m).u)
        fprintf(fid, '\t''%s''\n', ar.model(m).u{j});
    end
    fprintf(fid, '\t\t};\n');
    fprintf(fid,'\t\tu_cell{jc} = cell(size(u_labels));\n');

    
    
    %% switch over conditoins
    fprintf(fid,'\tswitch jc\n');

    for c = 1:length(ar.model(m).condition)
        xda = ~isempty(ar.model(m).x);
        uda = ~isempty(ar.model(m).u);
        
        c_p = ar.model(m).condition(c).pold;
        c_fp = ar.model(m).condition(c).fp;
        c_fp = ReplaceVarnames(c_fp);  % substitute each variable, e.g. x by p_struct.x
        
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
        fprintf(fid, '\t\t\tpcond_cell{jc} = {');
        for j=1:length(p)
            fprintf(fid, ' %s', p{j});
            if rem(j,10)==0
                fprintf(fid,'...\n\t\t\t\t');
            else
                fprintf(fid,' ');
            end
        end
        fprintf(fid, '};\n\n');

        fprintf(fid, '\t\t\tpcond_labels{jc} = {');
        for j=1:length(p)
            fprintf(fid, ' ''%s''', p{j});
            if rem(j,10)==0
                fprintf(fid,'...\n\t\t\t\t');
            else
                fprintf(fid,' ');
            end
        end
        fprintf(fid, '};\n\n');
       
        
        fprintf(fid,'\t\tx_struct = struct;\n');
        fprintf(fid,'\t\tx_cell{jc} = cell(size(x_labels));\n');
        if xda
%             fprintf(fid, '\t\t[t,x] = ode15s(@(t,x) %s(t,x,p,x0), [%f %f], x0);\n\n', fname, ar.model(m).tLim(1), ar.model(m).tLim(2));
            fprintf(fid, '\n\t\t[pred_cond{jc},x] = ode15s(@(t,x) %s(t,x,[pcond_cell{jc}{:}],x0), pred_cell{%i}, x0);\n', fname,ar.model(m).condition(c).dLink(1));
            fprintf(fid,'\t\tfor ix=1:length(x_labels)\n\t\t\tx_struct.(x_labels{ix}) = x(:,ix);\n\t\tend\n');
            fprintf(fid,'\t\tfor ix=1:length(x_labels)\n\t\t\tx_cell{jc}{ix} = x_struct.(x_labels{ix});\n\t\tend\n');
                
        else
            fprintf(fid, '\t\tpred_cond{jc} = pred_cell{%i};\n',ar.model(m).condition(c).dLink(1));
            fprintf(fid, '\t\tx=[];\n');
        end
            if uda            
            fprintf(fid, '\n');
            fprintf(fid,'\t\targsU = [fieldnames(p_struct)'',pcond_labels{jc},pred_labels];\n\n');
            fprintf(fid, '\t\tu = [');
            for j=1:length(u)
                fprintf(fid,'feval(inline(''%s'',argsU{:}),p_cell{:},pcond_cell{jc}{:},pred_cell{:});  ',strrep(strrep(strrep(u{j},'^','.^'),'/','./'),'*','.*'));%,pred_labels{ar.model(m).condition(c).dLink(1)});
            end
            fprintf(fid, '];\n');
            fprintf(fid,'\t\tif length(u)==1,\n');
            fprintf(fid,'\t\t\tu=u*ones(size(pred_cell{1}));\n');
            fprintf(fid,'\t\tend\n');

            fprintf(fid,'\t\tfor iu=1:length(u_labels)\n');
            fprintf(fid,'\t\t\tu_cond.(u_labels{iu}) = u(:,iu);\n');
            fprintf(fid,'\t\tend\n');

            fprintf(fid,'\t\tfor iu=1:length(u_labels)\n');
            fprintf(fid,'\t\t\tu_cell{jc}{iu} = u_cond.(u_labels{iu});\n');
            fprintf(fid,'\t\tend\n');
        else
            fprintf(fid, '\t\tu=[];\n');            
        end
    end
    
    fprintf(fid,'\n\tend %s switch over conditions\n','%');
    
    fprintf(fid,'\tz_struct=struct;\n');
    if(~isempty(ar.model(m).z)>0)
        fprintf(fid,'%%%% Derived variables:\n');
        fprintf(fid,'\targs{jc} = [fieldnames(x_struct)'',fieldnames(p_struct)'',pcond_labels{jc}];\n');
        for i=1:length(ar.model(m).z)            
            fprintf(fid,'\tz_struct.%s = feval(inline(''%s'',args{jc}{:}),x_cell{jc}{:},p_cell{:},pcond_cell{jc}{:});\n',ar.model(m).z{i},ar.model(m).fz{i});
        end
        fprintf(fid,'\tz_labels = fieldnames(z_struct);\n');
        fprintf(fid,'\tz_cell = cell(size(z_labels));\n');
        fprintf(fid,'\tfor i=1:length(z_labels)\n\t\tz_cell{i} = z_struct.(z_labels{i});\t\nend\n');

        fprintf(fid,'\n');
    else
        fprintf(fid,'\tz_labels = fieldnames(z_struct);\n');
        fprintf(fid,'\tz_cell = cell(size(z_labels));\n');
    end
    fprintf(fid,'\targs2{jc} = [fieldnames(x_struct)'',fieldnames(p_struct)'',fieldnames(z_struct)'',pcond_labels{jc},u_labels];\n');
    
    fprintf(fid, '%s plotting \n','%%');
    fprintf(fid, '\tif(~isempty(x))\n');
    fprintf(fid, '\t\tfigure;\n');
    fprintf(fid, '\t\tplot(pred_cond{jc},x);\n');
    fprintf(fid, '\t\tlegend(strrep(x_labels, ''_'', ''\\_''));\n');
    fprintf(fid, '\tend\n\n');

    fprintf(fid, '\tif(~isempty(u))\n');    
    fprintf(fid, '\t\tif length(u)==1,\n\t\t\t u=u*ones(size(pred_cell{%i}));\n\t\tend\n',ar.model(m).condition(c).dLink(1));
    fprintf(fid, '\t\tfigure;\n');
    fprintf(fid, '\t\tplot(pred_cell{%i},u);\n',ar.model(m).condition(c).dLink(1));
    fprintf(fid, '\t\tlegend(strrep(u_labels, ''_'', ''\\_''));\n');
    fprintf(fid, '\t\txlabel(strrep(pred_labels{jc},''_'',''\\_''));\n');    
    fprintf(fid, '\tend\n\n');
    
        
    if strcmp(version('-release'), '2016a')
        warning('on', 'symbolic:sym:sym:DeprecateExpressions')
    end
    
    %% save data as *.mat
    for d=1:length(ar.model(m).data)
        dat = struct;
        dat.tExp = ar.model(m).data(d).tExp;
        dat.tFine = ar.model(m).data(d).tFine;
        dat.yExp = ar.model(m).data(d).yExp;
        dat.yExpStd = ar.model(m).data(d).yExpStd;
        dat.fy   = ar.model(m).data(d).fy;
        dat.fystd   = ar.model(m).data(d).fystd;
        dat.dolog   = ar.model(m).data(d).logplotting;
        dat.m = m;
        dat.d = d;
        dat.cLink = ar.model(m).data(d).cLink;
        datfile = sprintf('%s_data_d%i.mat',ar.model(m).name,d);
        save(datfile,'dat');
        
        fprintf(fid,'\tload(''%s'',''dat'')\n',datfile);
        if d==1
            fprintf(fid,'\tdatasets = dat;\n',d);
        else
            fprintf(fid,'\tdatasets(%i) = dat;\n',d);
        end
        
        
    end
    fprintf(fid,'\tfor d=1:length(datasets);\n');
    fprintf(fid,'\t\tfor iy = 1:length(datasets(d).fy);\n');
    fprintf(fid,'\t\t\ty = feval(inline(datasets(d).fy{iy},args2{datasets(d).cLink}{:}),x_cell{datasets(d).cLink}{:},p_cell{:},z_cell{:},pcond_cell{jc}{:},u_cell{jc}{:});\n');
    fprintf(fid,'\t\t\tystd = feval(inline(datasets(d).fystd{iy},args2{datasets(d).cLink}{:}),x_cell{datasets(d).cLink}{:},p_cell{:},z_cell{:},pcond_cell{jc}{:},u_cell{jc}{:});\n');
    fprintf(fid,'\t\t\tfigure\n');
    fprintf(fid,'\t\t\tif datasets(d).dolog(iy)\n');
    fprintf(fid,'\t\t\t\tplot(pred_cond{datasets(d).cLink},log10(y),''k''), hold on\n');
    fprintf(fid,'\t\t\t\tplot(pred_cond{datasets(d).cLink},log10(y)+log10(ystd),''k--'')\n');
    fprintf(fid,'\t\t\t\tplot(pred_cond{datasets(d).cLink},log10(y)-log10(ystd),''k--'')\n');
    fprintf(fid,'\t\t\telse\n');
    fprintf(fid,'\t\t\t\tplot(pred_cond{datasets(d).cLink},y,''k''), hold on\n');
    fprintf(fid,'\t\t\t\tplot(pred_cond{datasets(d).cLink},y+ystd,''k--'')\n');
    fprintf(fid,'\t\t\t\tplot(pred_cond{datasets(d).cLink},y-ystd,''k--'')\n');
    fprintf(fid,'\t\t\tend\n');
    fprintf(fid,'\t\t\tif sum(~isnan(datasets(d).yExpStd(:,iy)))==0\n');
    fprintf(fid,'\t\t\t\tplot(datasets(d).tExp,datasets(d).yExp(:,iy),''o'')\n');
    fprintf(fid,'\t\t\telse\n');
    fprintf(fid,'\t\t\terrorbar(datasets(d).tExp,datasets(d).yExp,datasets(d).yExpStd,''o'')\n');
    fprintf(fid,'\t\t\tend\n');
    fprintf(fid,'\t\tend;\n');
    fprintf(fid,'\tend;\n');

    fprintf(fid,'end  % loop over sub-models\n');        
    fclose(fid);
    
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
        fp = strrep(char(fp),'pppp_','p_struct.');
    else
        fp = char(fp);
    end
end

function compile_arInputFunctionsC

