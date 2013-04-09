% Load Model definition
% 
% arLoadModel(name, m)
%
% name      filename of model definition file
% m         target position for model               [next]
%
% Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)

function arLoadModel(name, m)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

% load model from mat-file
if(~exist('Models','dir'))
    error('folder Models/ does not exist')
end
if(~exist(['Models/' name '.def'],'file'))
    error('model definition file %s.def does not exist in folder Models/', name)
end

if(~exist('m','var'))
    if(isfield(ar, 'model'))
        m = length(ar.model) + 1;
    else
        m = 1;
    end
end

fprintf('loading model #%i, from file Models/%s.def ...\n', m, name);
fid = fopen(['Models/' name '.def'], 'r');

% initial setup
ar.model(m).name = name;

% DESCRIPTION
str = textscan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(~strcmp(str{1},'DESCRIPTION'))
    error('parsing model %s for DESCRIPTION', name);
end
str = textscan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
ar.model(m).description = {};
while(~strcmp(str{1},'PREDICTOR'))
    ar.model(m).description(end+1,1) = str{1};
    str = textscan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
end

% PREDICTOR
C = textscan(fid, '%s %q %q %q %n %n\n',1, 'CommentStyle', ar.config.comment_string);
ar.model(m).t = cell2mat(C{1});
ar.model(m).tUnits(1) = C{2};
ar.model(m).tUnits(2) = C{3};
ar.model(m).tUnits(3) = C{4};
ar.model(m).tLim = [C{5} C{6}];
if(isnan(ar.model(m).tLim(1)))
    ar.model(m).tLim(1) = 0;
end
if(isnan(ar.model(m).tLim(2)))
    ar.model(m).tLim(2) = 10;
end

% COMPARTMENTS
ar.model(m).c = {};
ar.model(m).cUnits = {};
ar.model(m).pc = {};
ar.model(m).px = {};
C = textscan(fid, '%s %q %q %q %f\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'STATES'))
    if(~strcmp(C{1},'COMPARTMENTS'))
        ar.model(m).c(end+1) = C{1};
        ar.model(m).cUnits(end+1,1) = C{2};
        ar.model(m).cUnits(end,2) = C{3};
        ar.model(m).cUnits(end,3) = C{4};
        if(isnan(C{5}))
            ar.model(m).px(end+1) = {['vol_' cell2mat(C{1})]};
            ar.model(m).pc(end+1) = {['vol_' cell2mat(C{1})]};
        else
            ar.model(m).pc(end+1) = {num2str(C{5})};
        end
    end
    C = textscan(fid, '%s %q %q %q %f\n',1, 'CommentStyle', ar.config.comment_string);
end

% STATES
ar.model(m).px0 = {};
ar.model(m).x = {};
ar.model(m).xNames = {};
ar.model(m).xUnits = {};
ar.model(m).cLink = [];
ar.model(m).qPlotX = [];
C = textscan(fid, '%s %q %q %q %s %n %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'INPUTS'))
    if(length(cell2mat(C{1}))<2)
        error('STATE names need to be longer than 1');
    end
    ar.model(m).x(end+1) = C{1};
    ar.model(m).xUnits(end+1,1) = C{2};
    ar.model(m).xUnits(end,2) = C{3};
    ar.model(m).xUnits(end,3) = C{4};
    if(~isempty(ar.model(m).c))
        qcomp = ismember(ar.model(m).c, C{5});
        if(sum(qcomp)~=1)
            error('unknown compartement %s', cell2mat(C{5}));
        end
        ar.model(m).cLink(end+1) = find(qcomp);
    end
    if(isempty(C{6}) || isnan(C{6}))
        ar.model(m).qPlotX(end+1) = 1;
    else
        ar.model(m).qPlotX(end+1) = C{6};
    end
    if(~isempty(cell2mat(C{7})))
        ar.model(m).xNames(end+1) = C{7};
    else
        ar.model(m).xNames{end+1} = ar.model(m).x{end};
    end
    ar.model(m).px0(end+1) = {['init_' cell2mat(C{1})]};
    C = textscan(fid, '%s %q %q %q %s %n %q\n',1, 'CommentStyle', ar.config.comment_string);
end

% INPUTS
ar.model(m).u = {};
ar.model(m).uUnits = {};
ar.model(m).fu = {};
C = textscan(fid, '%s %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'REACTIONS') && ~strcmp(C{1},'REACTIONS-AMOUNTBASED') && ~strcmp(C{1},'ODES'))
    if(~strcmp(C{1},''))
        if(sum(ismember(ar.model(m).x, C{1}))>0)
            error('input %s already defined in STATES', cell2mat(C{1}));
        end
        ar.model(m).u(end+1) = C{1};
        ar.model(m).uUnits(end+1,1) = C{2};
        ar.model(m).uUnits(end,2) = C{3};
        ar.model(m).uUnits(end,3) = C{4};
        ar.model(m).fu(end+1,1) = C{5};
    end
    C = textscan(fid, '%s %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
end
ar.model(m).qPlotU = ones(size(ar.model(m).u));

% input parameters
varlist = cellfun(@symvar, ar.model(m).fu, 'UniformOutput', false);
ar.model(m).pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''});

% REACTIONS (respectively ODES)
ar.model(m).N = [];
ar.model(m).v = {};
ar.model(m).fv = {};
ar.model(m).fv_ma_reverse_pbasename = {};
ar.model(m).vUnits = {};
if(strcmp(C{1},'REACTIONS') || strcmp(C{1},'REACTIONS-AMOUNTBASED'))
    ar.model(m).isReactionBased = true;
    if(strcmp(C{1},'REACTIONS-AMOUNTBASED'))
        ar.model(m).isAmountBased = true;
    else
        ar.model(m).isAmountBased = false;
    end
    vcount = 1;
    str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
    while(~strcmp(str{1},'INVARIANTS'))
        source = {};
        while(~strcmp(str{1},'->') && ~strcmp(str{1},'<->'))
            if(~strcmp(str{1},'0') && ~strcmp(str{1},'+'))
                source(end+1) = str{1}; %#ok<AGROW>
            end
            str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
        end
        if(sum(~ismember(source, ar.model(m).x)) > 0)
            error('undefined source species in reaction %i: %s', vcount, ...
                source{~ismember(source, ar.model(m).x)})
        end
        
        if(strcmp(str{1},'<->'))
            reversible = true;
        else
            reversible = false;
        end
        
        target = {};
        str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
        while(~strcmp(str{1},'CUSTOM') && ~strcmp(str{1},'MASSACTION'))
            if(~strcmp(str{1},'0') && ~strcmp(str{1},'+'))
                target(end+1) = str{1}; %#ok<AGROW>
            end
            str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
        end
        if(sum(~ismember(target, ar.model(m).x)) > 0)
            error('undefined target species in reaction %i: %s', vcount, ...
                target{~ismember(target, ar.model(m).x)})
        end
        
        % infer flux units
        if(~isempty(source))
            ix = find(ismember(ar.model(m).x, source{1}));
        elseif(~isempty(target))
            ix = find(ismember(ar.model(m).x, target{1}));
        else
            error('reaction with empty N');
        end
        ar.model(m).vUnits{end+1,1} = [ar.model(m).xUnits{ix,1} '/' ar.model(m).tUnits{1}];
        ar.model(m).vUnits{end,2} = [ar.model(m).xUnits{ix,2} '/' ar.model(m).tUnits{2}];
        ar.model(m).vUnits{end,3} = [ar.model(m).xUnits{ix,3} '/' ar.model(m).tUnits{3}];
        
        if(strcmp(str{1},'MASSACTION'))
            massaction = true;
        else
            massaction = false;
        end
        
%         str = textscan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
        C = textscan(fid, '%q %q\n', 1, 'CommentStyle', ar.config.comment_string);
        str = C(1);
        ar.model(m).v{end+1} = cell2mat(C{2});
        
        ar.model(m).fv_ma_reverse_pbasename{end+1} = '';
        if(~massaction)
            if(reversible)
                error('reversible reactions for type CUSTOM not allowed.');
            else
                ar.model(m).fv(end+1,1) = str{1};
            end
        else
            if(~reversible)
                ar.model(m).fv{end+1,1} = cell2mat(str{1});
                for j=1:length(source)
                    ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' source{j}];
                end
            else
                ar.model(m).fv{end+1,1} = [cell2mat(str{1}) '_1'];
                ar.model(m).fv_ma_reverse_pbasename{end} = cell2mat(str{1});
                for j=1:length(source)
                    ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' source{j}];
                end
            end
        end
        
        % setup N
        ar.model(m).N(1:length(ar.model(m).x),vcount) = 0;
        for jj=1:length(source)
            for j=find(ismember(ar.model(m).x, source{jj}))
                ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) - 1;
            end
        end
        for jj=1:length(target)
            for j=find(ismember(ar.model(m).x, target{jj}))
                ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) + 1;
            end
        end
        
        % check for inconsistent educt compartments
        if(~isempty(ar.model(m).c) && ~ar.model(m).isAmountBased)
            for j=1:size(ar.model(m).N,2)
                if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)>0)))>1)
                    error('efflux from different compartments in reaction %s', ...
                        ar.model(m).fv{end});
                end
                if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)<0)))>1)
                    error('influx from different compartments in reaction %s', ...
                        ar.model(m).fv{end});
                end
            end
        end
        
        vcount = vcount + 1;
        
        % setup reversed reaction
        if(massaction && reversible)
            ar.model(m).fv{end+1,1} = [cell2mat(str{1}) '_2'];
            ar.model(m).fv_ma_reverse_pbasename{end+1} = cell2mat(str{1});
            for j=1:length(target)
                ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' target{j}];
            end
            
            % infer flux units
            if(~isempty(target))
                ix = find(ismember(ar.model(m).x, target{1}));
            elseif(~isempty(source))
                ix = find(ismember(ar.model(m).x, source{1}));
            else
                error('reaction with empty N');
            end
            ar.model(m).vUnits{end+1,1} = [ar.model(m).xUnits{ix,1} '/' ar.model(m).tUnits{1}];
            ar.model(m).vUnits{end,2} = [ar.model(m).xUnits{ix,2} '/' ar.model(m).tUnits{2}];
            ar.model(m).vUnits{end,3} = [ar.model(m).xUnits{ix,3} '/' ar.model(m).tUnits{3}];
            
            % setup N
            ar.model(m).N(1:length(ar.model(m).x),vcount) = 0;
            for jj=1:length(source)
                for j=find(ismember(ar.model(m).x, source{jj}))
                    ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) + 1;
                end
            end
            for jj=1:length(target)
                for j=find(ismember(ar.model(m).x, target{jj}))
                    ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) - 1;
                end
            end
            
            % check for inconsistent educt compartments
            if(~isempty(ar.model(m).c))
                for j=1:size(ar.model(m).N,2)
                    if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)>0)))>1)
                        error('efflux from different compartments in reaction %s', ...
                            ar.model(m).fv{end});
                    end
                    if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)<0)))>1)
                        error('influx from different compartments in reaction %s', ...
                            ar.model(m).fv{end});
                    end
                end
            end
            
            vcount = vcount + 1;
        end
        
        str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
    end
elseif(strcmp(C{1},'ODES'))
    ar.model(m).isReactionBased = false;
    str = textscan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
    ode_count = 0;
    while(~strcmp(str{1},'INVARIANTS'))
        if(~strcmp(str{1},''))
            ode_count = ode_count + 1;
            ar.model(m).fv{end+1,1} = cell2mat(str{1});
            ar.model(m).fv_ma_reverse_pbasename{end+1} = '';
            ar.model(m).vUnits{end+1,1} = [ar.model(m).xUnits{ode_count,1} '/' ar.model(m).tUnits{1}];
            ar.model(m).vUnits{end,2} = [ar.model(m).xUnits{ode_count,2} '/' ar.model(m).tUnits{2}];
            ar.model(m).vUnits{end,3} = [ar.model(m).xUnits{ode_count,3} '/' ar.model(m).tUnits{3}];
        end
        str = textscan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
    end
    if(ode_count ~= length(ar.model(m).x))
        error('number of ODES ~= number of variables');
    end
    ar.model(m).N = eye(length(ar.model(m).x));
end
ar.model(m).qPlotV = ones(1,length(ar.model(m).fv));
if(isempty(ar.model(m).fv))
    ar.model(m).isReactionBased = false;
end

% dynamic parameters
varlist = cellfun(@symvar, ar.model(m).fv, 'UniformOutput', false);
ar.model(m).px = union(setdiff(vertcat(varlist{:}), union(ar.model(m).x, ar.model(m).u)), ...
    ar.model(m).px);

% setup rhs
C = cell(size(ar.model(m).N));
if(length(ar.model(m).c)>1)    
    if(~isfield(ar.model(m),'isAmountBased') || ~ar.model(m).isAmountBased)
        for j=1:size(ar.model(m).N,1) % for every species j
            qinfluxwitheducts = ar.model(m).N(j,:) > 0 & sum(ar.model(m).N < 0,1) > 0;
            eductcompartment = zeros(size(qinfluxwitheducts));
            for jj=find(qinfluxwitheducts)
                eductcompartment(jj) = unique(ar.model(m).cLink(ar.model(m).N(:,jj)<0));
            end
            
            cfaktor = cell(size(qinfluxwitheducts));
            for jj=1:size(ar.model(m).N,2) % for every reaction jj
                if(qinfluxwitheducts(jj) && eductcompartment(jj)~=ar.model(m).cLink(j))
                    cfaktor{jj} = [ar.model(m).pc{eductcompartment(jj)} '/' ...
                        ar.model(m).pc{ar.model(m).cLink(j)}];
                else
                    cfaktor{jj} = '1';
                end
            end
            C(j,:) = transpose(cfaktor);
        end
    else
        for j=1:size(ar.model(m).N,1) % for every species j
            for jj=1:size(ar.model(m).N,2) % for every reaction jj
                C{j,jj} = ['1/' ar.model(m).pc{ar.model(m).cLink(j)}];
            end
        end
    end
else
    for j=1:size(ar.model(m).N,1) % for every species j
        for jj=1:size(ar.model(m).N,2) % for every reaction jj
            C{j,jj} = '1';
        end
    end
end
ar.model(m).fx = cell(length(ar.model(m).x),1);
tmpfx = (sym(ar.model(m).N).*sym(C)) * sym(ar.model(m).fv);
for j=1:length(ar.model(m).x) % for every species j
    ar.model(m).fx{j} = char(tmpfx(j));
end

% INVARIANTS
ar.model(m).fxeq = {};
C = textscan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'CONDITIONS'))
    if(~strcmp(C{1},''))
        ar.model(m).fxeq(end+1) = C{1};
    end
    C = textscan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
end

% extra invariational parameters
varlist = cellfun(@symvar, ar.model(m).fxeq, 'UniformOutput', false);
ar.model(m).pxeq = setdiff(vertcat(varlist{:}), union(ar.model(m).x, union(ar.model(m).u, ...
    ar.model(m).px)));

% collect parameters needed for ODE
ar.model(m).p = union(ar.model(m).px, union(ar.model(m).px0, union(ar.model(m).pc, union(ar.model(m).pxeq, ar.model(m).pu))));

% CONDITIONS
C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16);
ar.model(m).fp = transpose(ar.model(m).p);
while(~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS'))
    qcondpara = ismember(ar.model(m).p, C{1});
    if(sum(qcondpara)>0)
        ar.model(m).fp{qcondpara} = ['(' cell2mat(C{2}) ')'];
    else
        warning('unknown parameter in conditions: %s', cell2mat(C{1})); %#ok<WNTAG>
    end
    C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, ...
        'BufSize', 2^16-1);
end

% extra conditional parameters
varlist = cellfun(@symvar, ar.model(m).fp, 'UniformOutput', false);
ar.model(m).pcond = setdiff(setdiff(setdiff(vertcat(varlist{:}), ar.model(m).p), ar.model(m).x), ar.model(m).u);

% PARAMETERS
if(~isfield(ar, 'pExternLabels'))
    ar.pExternLabels = {};
    ar.pExtern = [];
    ar.qFitExtern = [];
    ar.qLog10Extern = [];
    ar.lbExtern = [];
    ar.ubExtern = [];
end
C = textscan(fid, '%s %f %n %n %f %f\n',1, 'CommentStyle', ar.config.comment_string);
while(~isempty(C{1}))
    ar.pExternLabels(end+1) = C{1};
    ar.pExtern(end+1) = C{2};
    ar.qFitExtern(end+1) = C{3};
    ar.qLog10Extern(end+1) = C{4};
    ar.lbExtern(end+1) = C{5};
    ar.ubExtern(end+1) = C{6};
    C = textscan(fid, '%s %f %n %n %f %f\n',1, 'CommentStyle', ar.config.comment_string);
end

fclose(fid);
