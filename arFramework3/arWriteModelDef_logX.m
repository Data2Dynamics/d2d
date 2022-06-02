%       zeroVal     which negative value should be used instead of log(0) ?
%                   [-30]
%                   either string or number
% 
% Example:
% modelname = 'FullModel';
% dataname = 'TimeCourseData';
% conv = arWriteModelDef_logX(ar.model,['Models/',modelname,'_logX.def']);
% arWriteDataDef_logX(ar.model.data,['Data/',dataname,'_logX.def'],conv);

function conv = arWriteModelDef_logX(model,file,zeroVal)
if ~exist('zeroVal','var') || isempty(zeroVal)
    zeroVal = '-30';
elseif isnumeric(zeroVal)
    zeroVal = num2str(zeroVal);
end
if(~exist('file','var') || isempty(file))
    file = [model.name,'_logX.def'];
end

states = model.x;
statesLG = strcat(model.x,'LG');
StatesLG = sym(statesLG);
statesLGexp = strcat('exp(',model.x,'LG)');
StatesLGexp = str2sym(statesLGexp);
States = sym(model.x);

statesInit = strcat('init_',model.x);
StatesInit = sym(statesInit);
StatesInitLGexp = str2sym(strcat('exp(init_',model.x,'LG)'));

conv.states = states;
conv.statesLG = statesLG;
conv.StatesLG = StatesLG;
conv.statesLGexp = statesLGexp;
conv.StatesLGexp = StatesLGexp
conv.States = States;
conv.statesInit = statesInit;
conv.StatesInit = StatesInit;
conv.StatesInitLGexp = StatesInitLGexp;


out = cell(0);
out{1} = 'DESCRIPTION';
for i=1:length(model.description)
    out{end+1} = sprintf('"%s"',model.description{i});
end
out{end+1} = sprintf('"arWriteModelDef_LogX.m: Automatic conversion to logX %s"\n',datestr(now));

out{end+1} = '';
out{end+1} = 'PREDICTOR';
out{end+1} = sprintf('%s\t%s\t"%s"\t"%s"\t%s\t%s',model.t,model.tUnits{1,1},model.tUnits{1,2},model.tUnits{1,3},num2str(model.tLim(1)),num2str(model.tLim(2)));
out{end+1} = '';
out{end+1} = 'COMPARTMENTS';
for i=1:length(model.c)
    out{end+1} = sprintf('%s\t%s\t"%s"\t"%s"\t%s',model.c{i},model.cUnits{i,1},model.cUnits{i,2},model.cUnits{i,3},model.pc{i});
end
out{end+1} = '';

out{end+1} = 'STATES';
if ~isempty(model.c)
    for i=1:length(model.x)
        % LG:
        out{end+1} = sprintf('%s\t%s[log]\t"%s[log]"\t"%s[log]"\t%s\t%i\t"%s"\t%s',statesLG{i},model.xUnits{i,1},model.xUnits{i,2},model.xUnits{i,3},model.c{model.cLink(i)}, model.qPlotX(i),statesLG{i}, '0');
%         out{end+1} = sprintf('%s\t%s\t"%s"\t"%s"\t%s\t%s\t"%s"\t%s',model.x{i},model.xUnits{i,1},model.xUnits{i,2},model.xUnits{i,3},model.c{1}, model.pc{1},model.xNames{i}, model.pc{model.c(i)});
    end
else
    for i=1:length(model.x)
        % LG:
        out{end+1} = sprintf('%s\t%s[log]\t"%s[log]"\t"%s[log]"\t%s\t%i\t"%s"\t%s',...
            statesLG{i},model.xUnits{i,1},model.xUnits{i,2},model.xUnits{i,3},'NONE', model.qPlotX(i),statesLG{i}, '0');
%         out{end+1} = sprintf('%s\t%sLG\t"%sLG"\t"%sLG"',statesLG{i},model.xUnits{i,1},model.xUnits{i,2},model.xUnits{i,3});
%         out{end+1} = sprintf('%s\t%s\t"%s"\t"%s"',model.x{i},model.xUnits{i,1},model.xUnits{i,2},model.xUnits{i,3});
    end
end

out{end+1} = '';
out{end+1}='INPUTS';
for i=1:length(model.fu)
    out{end+1} = sprintf('%s\t%s\t"%s"\t"%s"\t"%s"',...
        char(arSubs(arSubs(arSym(model.u{i}),States,StatesLGexp),StatesInit,StatesInitLGexp)),...
        model.uUnits{i,1},model.uUnits{i,2},model.uUnits{i,3},model.fu{i});
end
out{end+1} = '';

out{end+1}='REACTIONS';
%% ODEs
% for i=1:length(model.fx)
%     out{end+1} = sprintf('\t->\t%s\tCUSTOM\t"%s"',states{i},model.fx{i});
% end

% split fluxes: only ONE source or ONE target should occur for logX:
fluxes = cell(0);
fluxVs = cell(0);
source = cell(0);
target = cell(0);
for i=1:length(model.fv)
    if ~isempty(model.fv_source{i}) && ~isempty(model.fv_target{i})
        for j=1:length(model.fv_source{i})
            fluxes{end+1} = [model.Cm{j,i},'*',model.fv{i}];
            fluxVs{end+1} = ['v_',num2str(length(fluxes))];%model.v{i};
            source{end+1} = model.fv_source{i}(j);
            target{end+1} = cell(0);
        end
        
        for j=1:length(model.fv_target{i})
            fluxes{end+1} = [model.Cm{j,i},'*',model.fv{i}];
            fluxVs{end+1} = ['v_',num2str(length(fluxes))];%model.v{i};
            source{end+1} = cell(0);
            target{end+1} = model.fv_target{i}(j);
        end
        
    elseif ~isempty(model.fv_source{i})
        for j=1:length(model.fv_source{i})
            fluxes{end+1} = [model.Cm{j,i},'*',model.fv{i}];
            fluxVs{end+1} = ['v_',num2str(length(fluxes))];%model.v{i};
            source{end+1} = model.fv_source{i}(j);
            target{end+1} = cell(0);
        end
    elseif ~isempty(model.fv_target{i})
        for j=1:length(model.fv_target{i})
            fluxes{end+1} = [model.Cm{j,i},'*',model.fv{i}];
            fluxVs{end+1} = ['v_',num2str(length(fluxes))];%model.v{i};
            source{end+1} = cell(0);
            target{end+1} = model.fv_target{i}(j);
        end
    else
        error('case should not occur')
    end        
end

for i=1:length(fluxes)
    if ~isempty(source{i})
        sourceLG = arSubs(arSym(source{i}{:}),States,StatesLG);
        Fv_source = char(sourceLG);
        Fv_target = '';
        ind = strmatch(source{i},states,'exact');
        Fv = char(combine(arSubs(arSym(['(',fluxes{i},')/',char(statesLGexp(ind))]),States,StatesLGexp),'exp'));
    elseif ~isempty(target{i})
        targetLG = arSubs(arSym(target{i}{:}),States,StatesLG);
        Fv_source = '';
        Fv_target = char(targetLG);
        ind = strmatch(target{i},states,'exact');
        Fv = char(combine(arSubs(arSym(['(',fluxes{i},')/',char(statesLGexp(ind))]),States,StatesLGexp),'exp'));
    else
        error('Case should not occur');
    end
    
%     V = char(arSubs(arSym(fluxVs{i}),States,StatesLG));
    
    out{end+1} = sprintf('%s->%s\tCUSTOM\t"%s"\t"%s"',...
        sprintf('%s\t',Fv_source),...
        sprintf('\t%s',Fv_target),...
        Fv,...
        fluxVs{i});
end

out{end+1} = '';
out{end+1} = 'DERIVED';
%% Unlog states as derived
for i=1:length(model.x)
    out{end+1} = sprintf('%s\t%s\t%s\t%s\t%s\n',model.x{i}, model.xUnits{i,1}, model.xUnits{i,2},model.xUnits{i,3}, ['"exp(',statesLG{i},')"']);
end
%% Existing derived:
for i=1:length(model.z)
    out{end+1} = sprintf('%s\t%s\t%s\t%s\t"%s"',model.z{i},model.zUnits{i,1},model.zUnits{i,2},model.zUnits{i,3},model.fz{i});
end
out{end+1} = '';

%%
out{end+1} = 'OBSERVABLES';
for i=1:length(model.y)
    out{end+1} = sprintf('%s\t%s\t%s\t%s\t%i\t%i\t"%s"',model.y{i},model.yUnits{i,1},model.yUnits{i,2},model.yUnits{i,3},model.normalize(i),model.logfitting(i),model.fy{i});
end
out{end+1} = '';

out{end+1} = 'ERRORS';
for i=1:length(model.y)
    out{end+1} = sprintf('%s\t"%s"',model.y{i},model.fystd{i});
end
out{end+1} = '';

out{end+1} = 'CONDITIONS';
for i=1:length(model.fp)
    if strcmp(model.p{i},model.fp{i})~=1
        if ~isempty(strmatch(model.p{i},statesInit,'exact'))
            if strcmp(model.fp{i},'(0)')==1
                out{end+1} = sprintf('%s\t"%s"',[model.p{i},'LG'],zeroVal);
            else
                out{end+1} = sprintf('%s\t"%s"',[model.p{i},'LG'],...
                    ['log(',char(arSubs(arSubs(arSym(model.fp{i}),States,StatesLGexp),StatesInit,StatesInitLGexp)),')']);
            end
        else
            out{end+1} = sprintf('%s\t"%s"',...
                model.p{i},...
                char(arSubs(arSubs(arSym(model.fp{i}),States,StatesLGexp),StatesInit,StatesInitLGexp)));
        end
    end
end


fid = fopen(file,'w');
for i=1:length(out)
    fprintf(fid,'%s\n',out{i});
end
fclose(fid);

