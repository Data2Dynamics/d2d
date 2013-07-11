% export model to PottersWheel
%
% function arExportPW(m, c)
%
% m:    model index
% c:    condition index

function arExportPW(m, c)

global ar

if(~exist([cd '/PW' ], 'dir'))
    mkdir([cd '/PW' ])
end

fid = fopen(sprintf('./PW/%s.m', ar.model(m).name), 'w');

fprintf(fid, '%% PottersWheel model definition file\n\n');
fprintf(fid, 'function m = getModel() %%#ok<FNDEF>\n\n');

fprintf(fid, 'm             = pwGetEmptyModel();\n\n');

fprintf(fid, 'm.ID          = ''%s'';\n',  ar.model(m).name);
fprintf(fid, 'm.name        = ''%s'';\n',  ar.model(m).name);
fprintf(fid, 'm.description = ''%s'';\n',  ar.model(m).description{1});
fprintf(fid, 'm.authors     = {''%s''};\n',  ar.config.username);
fprintf(fid, 'm.dates       = {''%s''};\n',  datestr(now));
fprintf(fid, 'm.type        = ''PW-1-5'';\n\n');

% pwAddX
for jx = 1:length(ar.model(m).x)
    qp = ismember(ar.pLabel, ar.model(m).px0{jx}); %R2013a compatible
    if(sum(qp)==1)
        xvalue = ar.p(qp);
        xlb = ar.lb(qp);
        xub = ar.ub(qp);
        
        if(ar.qLog10(qp) == 1)
            xvalue = 10^xvalue;
            xlb = 10^xlb;
            xub = 10^xub;
        end
        
        if(ar.qFit(qp) == 1)
            typestr = '''global''';
        else
            typestr = '''fix''';
        end
            
        fprintf(fid, 'm = pwAddX(m, %20s, %20s, %10s, %20s, %20s, %20s, %20s, [], [], ''concentration'');\n',  ...
            sprintf('''%s''', ar.model(m).x{jx}), sprintf('%g', xvalue), typestr, sprintf('%g', xlb), sprintf('%g', xub), ...
            sprintf('''%s''', ar.model(m).xUnits{jx,2}), sprintf('''%s''', ar.model(m).c{ar.model(m).cLink(jx)}));
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).px0{jx}); %R2013a compatible
        xvalue = char(sym(ar.model(m).condition(c).fp{qp}));
        xlb = 0;
        xub = 1000;
        
        fprintf(fid, 'm = pwAddX(m, %20s, %20s, %10s, %20s, %20s, %20s, %20s, [], [], ''concentration'');\n',  ...
            sprintf('''%s''', ar.model(m).x{jx}), xvalue, '''fix''', sprintf('%g', xlb), sprintf('%g', xub), ...
            sprintf('''%s''', ar.model(m).xUnits{jx,2}), sprintf('''%s''', ar.model(m).c{ar.model(m).cLink(jx)}));
    else
        error('%s not found', ar.model(m).px0{jx});
    end
    

end
fprintf(fid, '\n');

% pwAddR
for jv = 1:length(ar.model(m).fv)
    sourcestr = '';
    count = 0;
    for jsource = find(ar.model(m).N(:,jv)<0)
        if(count>0)
            sourcestr = [sourcestr ', ']; %#ok<AGROW>
        end
        if(~isempty(jsource))
            sourcestr = [sourcestr '''' ar.model(m).x{jsource} '''']; %#ok<AGROW>
            count = count + 1;
        end
    end
    
    targetstr = '';
    count = 0;
    for jtarget = find(ar.model(m).N(:,jv)>0)
        if(count>0)
            targetstr = [targetstr ', ']; %#ok<AGROW>
        end
        if(~isempty(jtarget))
            targetstr = [targetstr '''' ar.model(m).x{jtarget} '''']; %#ok<AGROW>
            count = count + 1;
        end
    end
    
    ratetemplatestr = ar.model(m).fv{jv};
    ratetemplate = sym(ratetemplatestr);
    ratetemplate = subs(ratetemplate, ar.model(m).condition(c).pold, ...
        ar.model(m).condition(c).fp);
    
    count = 1;
    for jsource = find(ar.model(m).N(:,jv)<0)
        ratetemplate = subs(ratetemplate, sym(ar.model(m).x{jsource}), sym(sprintf('r%i', count)));
        count = count + 1;
    end
    count = 1;
    countp = 1;
    if(ratetemplate~=0)
        vars = symvar(ratetemplate);
        
        modifierstr = '';
        parameterstr = '';
        for jx = 1:length(vars)
            charvar = char(vars(jx));
            qx = ismember(ar.model(m).x, charvar); %R2013a compatible
            qu = ismember(ar.model(m).u, charvar); %R2013a compatible
            qp = ismember(ar.model(m).p, charvar); %R2013a compatible
            if(sum(qx)==1)
                if(count>1)
                    modifierstr = [modifierstr ', ']; %#ok<AGROW>
                end
                modifierstr = [modifierstr '''' ar.model(m).x{qx} '''']; %#ok<AGROW>
                ratetemplate = subs(ratetemplate, vars(jx), sym(sprintf('m%i', count)));
                count = count + 1;
            elseif(sum(qx)>1)
                error('multiple species matches error');
            elseif(sum(qu)==1)
                if(count>1)
                    modifierstr = [modifierstr ', ']; %#ok<AGROW>
                end
                modifierstr = [modifierstr '''' ar.model(m).u{qu} '''']; %#ok<AGROW>
                ratetemplate = subs(ratetemplate, vars(jx), sym(sprintf('m%i', count)));
                count = count + 1;
            elseif(sum(qx)>1)
                error('multiple species matches error');
            elseif(sum(qp)==1)
                if(countp>1)
                    parameterstr = [parameterstr ', ']; %#ok<AGROW>
                end
                parameterstr = [parameterstr '''' ar.model(m).p{qp} '''']; %#ok<AGROW>
                ratetemplate = subs(ratetemplate, vars(jx), sym(sprintf('k%i', countp)));
                countp = countp + 1;
            elseif(sum(qp)>1)
                error('multiple parameter matches error');
            end
        end
        
        
        fprintf(fid, 'm = pwAddR(m, {%s}, {%s}, {%s}, ''C'', [], %s, {%s});\n', ...
            sourcestr, targetstr, modifierstr, ['''' char(ratetemplate) ''''], parameterstr);
    end
end
fprintf(fid, '\n');

% pwAddC
for jc = 1:length(ar.model(m).c)
    qp = ismember(ar.pLabel, ar.model(m).pc{jc}); %R2013a compatible
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        fprintf(fid, 'm = pwAddC(m, %20s, %20s, [], 3, %20s, %20s, 1);\n',  ...
            sprintf('''%s''', ar.model(m).c{jc}), sprintf('%g', 10^pvalue), ...
            sprintf('''%s''', ar.model(m).c{jc}), sprintf('''%s''', ar.model(m).xUnits{jx,2}));
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).pc{jc}); %R2013a compatible
        if(sum(qp)==1)
            pvalue = ar.model(m).condition(c).fp{qp};
            fprintf(fid, 'm = pwAddC(m, %20s, %20s, [], 3, %20s, %20s, 1);\n',  ...
                sprintf('''%s''', ar.model(m).c{jc}), pvalue, ...
                sprintf('''%s''', ar.model(m).c{jc}), sprintf('''%s''', ar.model(m).xUnits{jx,2}));
        else
            error('%s not found', ar.model(m).pc{jc});
        end
    else
        error('%s not found', ar.model(m).pc{jc});
    end
    
    
end
fprintf(fid, '\n');

% pwAddK
for jp = 1:length(ar.model(m).condition(c).p)
    qp = ismember(ar.pLabel, ar.model(m).condition(c).p{jp}); %R2013a compatible
    
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        plb = ar.lb(qp);
        pub = ar.ub(qp);
        
        if(ar.qLog10(qp) == 1)
            pvalue = 10^pvalue;
            plb = 10^plb;
            pub = 10^pub;
        end
        
        if(ar.qFit(qp) == 1)
            typestr = '''global''';
        else
            typestr = '''fix''';
        end
        
    else
        pvalue = 1;
        plb = 0;
        pub = 1000;
        
        typestr = '''fix''';
    end
    
    fprintf(fid, 'm = pwAddK(m, %20s, %20s, %10s, %20s, %20s, [], [], []);\n',  ...
        sprintf('''%s''', ar.model(m).condition(c).p{jp}), sprintf('%g', pvalue), ...
        typestr, sprintf('%g', plb), sprintf('%g', pub));
end
fprintf(fid, '\n');

% pwAddU
for ju = 1:length(ar.model(m).u)
    qp = ismember(ar.pLabel, ar.model(m).condition(c).fu{ju}); %R2013a compatible
    
    if(sum(qp)==1)
        xvalue = ar.p(qp);
        fprintf(fid, 'm = pwAddU(m, %20s, ''steps'', [-100 0], [0 %s]);\n', ...
            sprintf('''%s''', ar.model(m).u{ju}), sprintf('%g', 10^xvalue));
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).condition(c).fu{ju}); %R2013a compatible
        xvalue = ar.model(m).condition(c).fp{qp};
        fprintf(fid, 'm = pwAddU(m, %20s, ''steps'', [-100 0], [0 %s]);\n', ...
            sprintf('''%s''', ar.model(m).u{ju}), xvalue);
    end
end
fprintf(fid, '\n');

fprintf(fid, 'm.t = %g:%g:%g;\n\n', ar.model(m).condition(c).tFine(1), ...
    ar.model(m).condition(c).tFine(2)-ar.model(m).condition(c).tFine(1), ar.model(m).condition(c).tFine(end));


fclose(fid);
