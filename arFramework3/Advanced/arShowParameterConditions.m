function arShowParameterConditions(m,d)

global ar

if(~exist('m','var'))
    m = 1;
end
if(~exist('d','var'))
    d = 1;
end

maxlabellength = max(cellfun(@length, ar.model(m).p));
fprintf('Model m%i, %s\n', m, ar.model(m).name);
for j=1:length(ar.model(m).p)
    if(~strcmp(ar.model(m).p{j},ar.model(m).fp{j}))
        fprintf('\t%s -> %s\n', arExtendStr(ar.model(m).p{j},maxlabellength), ar.model(m).fp{j});
    end
end
fprintf('\n');

if(isfield(ar.model(m), 'data'))
    if(isfield(ar.model(m).data(d),'pold'))
        maxlabellength = max(cellfun(@length, ar.model(m).data(d).pold));
        fprintf('Data d%i, %s\n', d, ar.model(m).data(d).name);
        for j=1:length(ar.model(m).data(d).pold)
            if(~strcmp(ar.model(m).data(d).pold{j},ar.model(m).data(d).fp{j}))
                fprintf('\t%s -> %s\n', arExtendStr(ar.model(m).data(d).pold{j},maxlabellength), ar.model(m).data(d).fp{j});
            end
        end
        fprintf('\n');
        
        c = ar.model(m).data(d).cLink;
        maxlabellength = max(cellfun(@length, ar.model(m).condition(c).p));
        fprintf('Condition c%i\n', c);
        for j=1:length(ar.model(m).condition(c).pold)
            if(~strcmp(ar.model(m).condition(c).pold{j},ar.model(m).condition(c).fp{j}))
                fprintf('\t%s -> %s\n', arExtendStr(ar.model(m).condition(c).pold{j},maxlabellength), ar.model(m).condition(c).fp{j});
            end
        end
        fprintf('\n');
    else
        maxlabellength = max(cellfun(@length, ar.model(m).data(d).p));
        fprintf('Data d%i, %s\n', d, ar.model(m).data(d).name);
        for j=1:length(ar.model(m).data(d).p)
            if(~strcmp(ar.model(m).data(d).p{j},ar.model(m).data(d).fp{j}))
                fprintf('\t%s -> %s\n', arExtendStr(ar.model(m).data(d).p{j},maxlabellength), ar.model(m).data(d).fp{j});
            end
        end
        fprintf('\n');
    end
else
    if(isfield(ar.model(m), 'condition'))
        if(isfield(ar.model(m).condition,'pold'))
            maxlabellength = max(cellfun(@length, ar.model(m).condition.pold));
            fprintf('Condition c%i\n', 1);
            for j=1:length(ar.model(m).condition.pold)
                if(~strcmp(ar.model(m).condition.pold{j},ar.model(m).condition.fp{j}))
                    fprintf('\t%s -> %s\n', arExtendStr(ar.model(m).condition.pold{j},maxlabellength), ar.model(m).condition.fp{j});
                end
            end
            fprintf('\n');
        else
            maxlabellength = max(cellfun(@length, ar.model(m).condition.p));
            fprintf('Condition c%i\n', 1);
            for j=1:length(ar.model(m).condition.p)
                if(~strcmp(ar.model(m).condition.p{j},ar.model(m).condition.fp{j}))
                    fprintf('\t%s -> %s\n', arExtendStr(ar.model(m).condition.p{j},maxlabellength), ar.model(m).condition.fp{j});
                end
            end
            fprintf('\n');
        end
    end
end