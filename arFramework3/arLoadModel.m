% Load Model definition to next free index position
% 
% arLoadModel(name)
%
% name      filename of model definition file
% m         model index for the model to go into
%
% 'ModelPath'            Path to the data files.
%                       Default: DataPath = 'Data/'
% 
% Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)

function arLoadModel(name, m, varargin)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

% Name for a state which is added to a stateless model to make sure that
% the required fields in the c-code aren't empty.
NOSTATE = '___dummy___';

% custom
switches = { 'conditions', 'modelpath' };
extraArgs = [ 1, 1 ];
description = { ...
    {'', 'Specified extra conditions'}, {'', 'Path to the model definition files'} };
    
opts = argSwitch( switches, extraArgs, description, 1, varargin );
if isempty(opts.modelpath_args)
    ModelPath = 'Models/';
else
    ModelPath = strtrim(opts.modelpath_args);
    if ModelPath(end)~='/' && ModelPath(end)~='\'
        ModelPath = [ModelPath,'/'];
    end
end


% load model from mat-file
if(~exist(ModelPath,'dir'))
    error('folder %s does not exist',ModelPath)
end
if strcmp(strrep(name,' ',''),name)~=1
    name
    error('File names should not contain empty spaces. Please remove it.');
end
if(~exist([ModelPath,  name '.def'],'file'))
    error('model definition file %s.def does not exist in folder %s', name, ModelPath)
end

if(~exist('m','var')) || isempty(m)
    if(isfield(ar, 'model'))
        m = length(ar.model) + 1;
    else
        m = 1;
    end
else
    error(['Usage arLoadModel(name, m) is deprecated. Please use arLoadModel(name) ' ...
        'and note that the model will be loaded to the next free index position by default.']);
end

% remember the function call
ar.setup.commands{end+1} = mfilename; % this file name
ar.setup.arguments{end+1} = {name, [],varargin{:}}; % argument m is deprecated and is therefore set as empty
ar.setup.datafiles{end+1} = {};
ar.setup.modelfiles{end+1} = [ModelPath, name,'.def'];

% Disable this if you are having problems because of the preprocessor
preprocessor = 1;
arFprintf(1, 'loading model #%i, from file %s%s.def...\n', m, ModelPath, name);
if ( ~preprocessor )
    fid = fopen([ModelPath,  name '.def'], 'r');
else
    % Load into a struct
    fid.fn  = [ModelPath,  name '.def'];
    fid.str = fileread([ModelPath,  name '.def']);
    fid.pos = 1;
    
    fid = arPreProcessor(fid);
end

% initial setup
ar.model(m).name = name;
ar.model(m).path = [pwd,filesep,ModelPath];

% Validate input
if ( opts.conditions )
    if ~(size( opts.conditions_args, 2 ) == 2) || ~iscell( opts.conditions_args )
        error( 'Additional conditions should be given in Mx2 cell array.' );
    else
        ar.model(m).extra_conditions = opts.conditions_args;
    end
else
    ar.model(m).extra_conditions = {};
end

matVer = arVer;

% DESCRIPTION
[str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(isempty(strfind(str{1},'DESCRIPTION')))
    arParsingError( fid, 'parsing model %s for DESCRIPTION', ar.model(m).name);
end

% check version
if(strcmp(str{1},'DESCRIPTION'))
    % def_version = 1;
elseif(strcmp(str{1},'DESCRIPTION-V2'))
    arParsingError( fid, 'DESCRIPTION-V2 not supported yet');
else
    arParsingError( fid, 'invalid version identifier: %s', cell2mat(str{1}));
end

% read comments
[str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
ar.model(m).description = {};
while(~strcmp(str{1},'PREDICTOR'))
    ar.model(m).description(end+1,1) = str{1};
    [str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
end

% PREDICTOR
[C, fid] = arTextScan(fid, '%s %q %q %q %n %n\n',1, 'CommentStyle', ar.config.comment_string);
arValidateInput( C, 'predictor', 'identifier for independent variable', 'unit type', 'unit', 'label for plotting' );
ar.model(m).t = cell2mat(C{1});
ar.model(m).tUnits(1) = C{2};
ar.model(m).tUnits(2) = C{3};
ar.model(m).tUnits(3) = C{4};
ar.model(m).tLim = [checkNum(C{5}, 0) checkNum(C{6}, 10)];

% COMPARTMENTS
ar.model(m).c = {};
ar.model(m).cUnits = {};
ar.model(m).pc = {};
ar.model(m).px = {};
[C, fid] = arTextScan(fid, '%s %q %q %q %f\n',1, 'CommentStyle', ar.config.comment_string);
if(~strcmp(C{1},'COMPARTMENTS'))
    arParsingError( fid, 'currently only one PREDICTOR allowed');
end
while(~strcmp(C{1},'STATES'))
    if(~strcmp(C{1},'COMPARTMENTS'))
        arValidateInput( C, 'compartments', 'compartment', 'unit type (i.e. V)', 'unit (i.e. pl)', 'label (i.e. "vol.")' );
        ar.model(m).c(end+1) = C{1};
        ar.model(m).cUnits(end+1,1) = C{2};
        ar.model(m).cUnits(end,2) = C{3};
        ar.model(m).cUnits(end,3) = C{4};
        
        if(isempty(C{5})||isnan(C{5}))
            ar.model(m).px(end+1) = {['vol_' cell2mat(C{1})]};
            ar.model(m).pc(end+1) = {['vol_' cell2mat(C{1})]};
        else
            ar.model(m).pc(end+1) = {num2str(C{5})};
        end
    end
    [C, fid] = arTextScan(fid, '%s %q %q %q %f\n',1, 'CommentStyle', ar.config.comment_string);
end

% STATES
ar.model(m).px0 = {};
ar.model(m).x = {};
ar.model(m).xNames = {};
ar.model(m).xUnits = {};
ar.model(m).cLink = [];
ar.model(m).qPlotX = [];
ar.model(m).qPositiveX = [];
[C, fid] = arTextScan(fid, '%s %q %q %q %s %n %q %n\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'INPUTS'))
    if ( strcmp( C{1}, 'REACTIONS' ) )
        arParsingError( fid,  'Missing field INPUTS. This section should be specified after STATES and before REACTIONS. See: "Setting up models"' );
    end    
    
    arValidateInput( C, 'state', 'unique identifier', 'unit type (i.e. C)', 'unit (i.e. nM)', 'label for plots (i.e. "conc.")' );
    
    if(length(cell2mat(C{1}))<2)
        arParsingError( fid, 'STATE names need to be longer than 1');
    end
    try
        if(isempty(symvar(arMyStr2Sym(C{1}))))
            arParsingError( fid, 'STATE name ''%s'' is reserved by MATLAB. Please rename!',cell2mat(C{1}));
        end
    catch
        C
        C{1}
        arParsingError( fid, 'Malformed input in STATE section.' );
    end
    
    ar.model(m).x(end+1) = C{1};
    ar.model(m).xUnits(end+1,1) = C{2};
    ar.model(m).xUnits(end,2) = C{3};
    ar.model(m).xUnits(end,3) = C{4};
    if(~isempty(ar.model(m).c))
        qcomp = ismember(ar.model(m).c, C{5}); %R2013a compatible
        if(sum(qcomp)~=1)
            arParsingError( fid, 'unknown compartment %s', cell2mat(C{5}));
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
    if(isempty(C{8}) || isnan(C{8}))
        ar.model(m).qPositiveX(end+1) = 0;
    else
        ar.model(m).qPositiveX(end+1) = C{8};
    end
    ar.model(m).px0(end+1) = {['init_' cell2mat(C{1})]};
    [C, fid] = arTextScan(fid, '%s %q %q %q %s %n %q %n\n',1, 'CommentStyle', ar.config.comment_string);
end

% Workaround for models that only use observation functions
if ( isempty( ar.model(m).x ) )
	ar.model(m).x(end+1) = {NOSTATE};
    ar.model(m).xUnits(end+1,1) = {'NA'};
    ar.model(m).xUnits(end,2) = {'NA'};
    ar.model(m).xUnits(end,3) = {'NA'};
    ar.model(m).cLink(end+1) = 0;
    ar.model(m).qPlotX(end+1) = 0;
    ar.model(m).xNames(end+1) = {NOSTATE};
    ar.model(m).qPositiveX(end+1) = 0;
    ar.model(m).px0(end+1) = {['init_' NOSTATE]};
    
    warning( 'D2D does not support models without state variables or equations. Adding species to model.' );
end

% INPUTS
ar.model(m).u = {};
ar.model(m).uUnits = {};
ar.model(m).fu = {};
ar.model(m).uNames = {};
[C, fid] = arTextScan(fid, '%s %q %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'REACTIONS') && ~strcmp(C{1},'REACTIONS-AMOUNTBASED') && ~strcmp(C{1},'ODES'))
    if(~strcmp(C{1},''))
        arValidateInput( C, 'input', 'unique input name', 'unit type (i.e. C)', 'unit (i.e. "units/cell")', 'plain text label for plots ("conc.")', 'input function' );
        if(sum(ismember(ar.model(m).x, C{1}))>0) %R2013a compatible
            arParsingError( fid, 'input %s already defined in STATES', cell2mat(C{1}));
        end
        ar.model(m).u(end+1) = C{1};
        ar.model(m).uUnits(end+1,1) = C{2};
        ar.model(m).uUnits(end,2) = C{3};
        ar.model(m).uUnits(end,3) = C{4};
        ar.model(m).fu(end+1,1) = C{5};
        if(~isempty(cell2mat(C{6})))
            ar.model(m).uNames(end+1) = C{6};
        else
            ar.model(m).uNames{end+1} = '';
        end
    end
    [C, fid] = arTextScan(fid, '%s %q %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
end
ar.model(m).qPlotU = ones(size(ar.model(m).u));

% input parameters
varlist = cellfun(@symvar, ar.model(m).fu, 'UniformOutput', false);
ar.model(m).pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}); %R2013a compatible

% REACTIONS (or ODES)
ar.model(m).N = [];
ar.model(m).v = {};
ar.model(m).fv = {};
ar.model(m).fv_source = {};
ar.model(m).fv_target = {};
ar.model(m).fv_sourceCoeffs = {};
ar.model(m).fv_targetCoeffs = {};
ar.model(m).reversible = [];
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
    %str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
    
    % Read single line (easier to trace errors back to their line)
    [ line, remainder, fid ] = readLine( fid, ar.config.comment_string );
    %str = textscan(line, '%s', 1);
    [str, remainder] = grabtoken( remainder, '%s', 1 );
    while(~strcmp(str{1},'INVARIANTS') && ~strcmp(str{1},'DERIVED'))
        source = {};
        sourceCoeffs = [];

        if ( strcmp(str{1}, 'OBSERVABLES') || strcmp(str{1}, 'CONDITIONS') )
            arParsingError( fid, 'Missing field DERIVED. This section should be specified after REACTIONS and before OBSERVABLES / CONDITIONS. See: "Setting up models paragraph 1.7"');
        end
        nextValue = 1;
        while(~strcmp(str{1},'->') && ~strcmp(str{1},'<->'))
            if(~strcmp(str{1},'0') && ~strcmp(str{1},'+'))
                % Check whether a stoichiometric coefficient is specified
                if ~isempty(str2num(str{1}{1})) %#ok
                    nextValue = str2num(str{1}{1}); %#ok
                else
                    source(end+1) = str{1}; %#ok<AGROW>
                    sourceCoeffs(end+1) = nextValue; %#ok<AGROW>
                    nextValue = 1;
                end
            end
            %str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
            [str, remainder] = grabtoken( remainder, '%s', 1 );
            if ( isempty(str{1}) )
                arParsingError( fid, 'incomplete reaction definition in reaction %i: %s', vcount, line)
            end
        end
        if(sum(~ismember(source, ar.model(m).x)) > 0) %R2013a compatible
            arParsingError( fid, '%s\nundefined source species in reaction %i: %s', line, vcount, ...
                source{~ismember(source, ar.model(m).x)}) %R2013a compatible
        end
        
        if(strcmp(str{1},'<->'))
            reversible = true;
        else
            reversible = false;
        end
        
        ar.model(m).reversible(end+1) = reversible;
        
        target = {};
        targetCoeffs = [];
        nextValue = 1;
        %str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
        [str, remainder] = grabtoken( remainder, '%s', 1 );
        while(~strcmp(str{1},'CUSTOM') && ~strcmp(str{1},'MASSACTION') && ~strcmp(str{1},'MASSACTIONKD'))
            if(~strcmp(str{1},'0') && ~strcmp(str{1},'+'))
                % Check whether a stoichiometric coefficient is specified
                if ( ~isempty(str2num(str{1}{1})) )  %#ok
                    nextValue = str2num(str{1}{1});   %#ok
                else
                    target(end+1)       = str{1}; %#ok<AGROW> 
                    targetCoeffs(end+1) = nextValue; %#ok<AGROW>
                    nextValue = 1;
                end
            end
            %str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
            [str, remainder] = grabtoken(remainder, '%s', 1);
            
            if ( isempty( str{1} ) )
                arParsingError( fid, 'missing keyword CUSTOM, MASSACTION or MASSACTIONKD before reaction rate expression');
            end
        end
        if(sum(~ismember(target, ar.model(m).x)) > 0) %R2013a compatible
            arParsingError( fid, 'undefined target species in reaction %i: %s', vcount, ...
                target{~ismember(target, ar.model(m).x)}) %R2013a compatible
        end
        
        % infer flux units
        if(~isempty(source))
            ix = find(ismember(ar.model(m).x, source{1})); %R2013a compatible
        elseif(~isempty(target))
            ix = find(ismember(ar.model(m).x, target{1})); %R2013a compatible
        else
            arParsingError( fid, 'reaction with empty N');
        end
        ar.model(m).vUnits{end+1,1} = [ar.model(m).xUnits{ix,1} '/' ar.model(m).tUnits{1}];
        ar.model(m).vUnits{end,2} = [ar.model(m).xUnits{ix,2} '/' ar.model(m).tUnits{2}];
        ar.model(m).vUnits{end,3} = [ar.model(m).xUnits{ix,3} '/' ar.model(m).tUnits{3}];
        
        if(strcmp(str{1},'MASSACTION'))
            massaction = true;
            massactionkd = false;
        else
            massaction = false;
        end
        if(strcmp(str{1},'MASSACTIONKD'))
            massaction = true;
            massactionkd = true;
        end
        
        [C, remainder] = grabtoken(remainder, '%q %q', 1);
        
        %C = arTextScan(fid, '%q %q\n', 1, 'CommentStyle', ar.config.comment_string);
        arValidateInput(C, 'REACTIONS', 'reaction rate expression' );
        str = C(1);
        if ( ~isempty(C{2}) )
            ar.model(m).v{end+1} = cell2mat(C{2});
        else
            ar.model(m).v{end+1} = sprintf('v_%d', length(ar.model(m).v) );
        end
        
        ar.model(m).fv_ma_reverse_pbasename{end+1} = '';
        if(~massaction)
            if(reversible)
                warning('Reversible reactions for type CUSTOM in experimental phase. Proceed with caution!');
            end
            ar.model(m).fv(end+1,1) = str{1};
            
            % check for negative fluxes possible
            if(ar.config.checkForNegFluxes==1)
                if (~reversible)
                    try
                        symtmp = arMyStr2Sym(str{1});
                    catch
                        if ( iscell( str{1} ) )
                            arParsingError( fid,  'Parsing error in REACTIONS at %s in model %d', str{1}{:}, m );
                        else
                            arParsingError( fid,  'Parsing error in REACTIONS at %s in model %d', str{1}, m );
                        end
                    end
                    for j=1:length(source)
                        symtmpsubs = subs(symtmp, arMyStr2Sym(source{j}), 0);
                        sva = symvar(symtmp);
                        symtmpsubs = subs(symtmpsubs, sva, rand(1, numel(sva)) );
                        if(symtmpsubs~=0)                           
                            arFprintf(1, 2, 'Possible negative flux in reaction #%i:\n', length(ar.model(m).fv));
                            arFprintf(1, 2, '%s : %s\n', arAssembleReactionStr(source, target, false, sourceCoeffs, targetCoeffs), cell2mat(str{1}));
                            arFprintf(1, 2, 'Source species %s missing ?\n\n', source{j});
                            arFprintf(1, 2, 'Deactivate this error message with: ar.config.checkForNegFluxes = false;\n\n');
                            arParsingError( fid, 'Possible negative fluxes in reaction');
                        end
                    end
                else
                    symtmp = arMyStr2Sym(str{1});
                    for j=1:length(source)
                        symtmpsubs = subs(symtmp, arMyStr2Sym(source{j}), 0);
                        for k=1:length(target)
                            symtmpsubs = subs(symtmpsubs, arMyStr2Sym(target{k}), 0);
                        end
                        sva = symvar(symtmp);
                        symtmpsubs = subs(symtmpsubs, sva, rand(1, numel(sva)) );
                        
                        if(symtmpsubs~=0)
                            arFprintf(1, 2, 'Possible flux in reaction without presence of source #%i:\n', length(ar.model(m).fv));
                            arFprintf(1, 2, '%s : %s\n', arAssembleReactionStr(source, target, false, sourceCoeffs, targetCoeffs), cell2mat(str{1}));
                            arFprintf(1, 2, 'Source species %s missing ?\n\n', source{j});
                            arFprintf(1, 2, 'Deactivate this error message with: ar.config.checkForNegFluxes = false;\n\n');
                            warning('Possible erroneous fluxes in reaction');
                        end                        
                    end
                    for j=1:length(target)
                        symtmpsubs = subs(symtmp, arMyStr2Sym(target{j}), 0);
                        for k=1:length(source)
                            symtmpsubs = subs(symtmpsubs, arMyStr2Sym(source{k}), 0);
                        end
                        
                        if(symtmpsubs~=0)
                            arFprintf(1, 2, 'Possible flux in reaction without presence of product #%i:\n', length(ar.model(m).fv));
                            arFprintf(1, 2, '%s : %s\n', arAssembleReactionStr(source, target, false, sourceCoeffs, targetCoeffs), cell2mat(str{1}));
                            arFprintf(1, 2, 'Product species %s missing ?\n\n', target{j});
                            arFprintf(1, 2, 'Deactivate this error message with: ar.config.checkForNegFluxes = false;\n\n');
                            warning('Possible erroneous fluxes in reaction');
                        end                        
                    end                    
                end
            end
        else
            if(~reversible)
                ar.model(m).fv{end+1,1} = cell2mat(str{1});
                for j=1:length(source)
                    if ( sourceCoeffs ~= 1 )
                        ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' source{j} '^' num2str(sourceCoeffs(j))];
                    else
                        ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' source{j} ];
                    end
                end
            else
                if(massactionkd)
                    ar.model(m).fv{end+1,1} = [cell2mat(str{1}) '_1*' cell2mat(str{1}) '_2'];
                else
                    ar.model(m).fv{end+1,1} = [cell2mat(str{1}) '_1'];
                end
                ar.model(m).fv_ma_reverse_pbasename{end} = cell2mat(str{1});
                for j=1:length(source)
                    if ( sourceCoeffs ~= 1 )
                        ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' source{j} '^' num2str(sourceCoeffs(j))];
                    else
                        ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' source{j} ];
                    end
                end
            end
        end
        
        % setup N
        ar.model(m).N(1:length(ar.model(m).x),vcount) = 0;
        for jj=1:length(source)
            for j=find(ismember(ar.model(m).x, source{jj})) %R2013a compatible
                ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) - sourceCoeffs(jj);
            end
        end
        for jj=1:length(target)
            for j=find(ismember(ar.model(m).x, target{jj})) %R2013a compatible
                ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) + targetCoeffs(jj);
            end
        end
        ar.model(m).fv_source{end+1,1} = source;
        ar.model(m).fv_target{end+1,1} = target;
        ar.model(m).fv_sourceCoeffs{end+1,1} = sourceCoeffs;
        ar.model(m).fv_targetCoeffs{end+1,1} = targetCoeffs;        
        
        % check for inconsistent educt compartments
        if(~isempty(ar.model(m).c) && ~ar.model(m).isAmountBased)
            for j=1:size(ar.model(m).N,2)
                if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)>0)))>1)
                    arParsingError( fid, 'efflux from different compartments in reaction %s', ...
                        ar.model(m).fv{end});
                end
                if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)<0)))>1)
                    arParsingError( fid, 'influx from different compartments in reaction %s', ...
                        ar.model(m).fv{end});
                end
            end
        end
        
        vcount = vcount + 1;
        
        % setup reversed reaction
        if(massaction && reversible)
            ar.model(m).fv{end+1,1} = [cell2mat(str{1}) '_2'];
            ar.model(m).fv_source{end+1,1} = ar.model(m).fv_target{end,1};
            ar.model(m).fv_target{end+1,1} = ar.model(m).fv_source{end,1};
            ar.model(m).fv_ma_reverse_pbasename{end+1} = cell2mat(str{1});
            for j=1:length(target)
                if ( targetCoeffs(j) ~= 1 )
                    ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' target{j} '^' num2str(targetCoeffs(j))];
                else
                    ar.model(m).fv{end,1} = [ar.model(m).fv{end,1} '*' target{j}];
                end
            end
            
            % infer flux units
            if(~isempty(target))
                ix = find(ismember(ar.model(m).x, target{1})); %R2013a compatible
            elseif(~isempty(source))
                ix = find(ismember(ar.model(m).x, source{1})); %R2013a compatible
            else
                arParsingError( fid, 'reaction with empty N');
            end
            ar.model(m).vUnits{end+1,1} = [ar.model(m).xUnits{ix,1} '/' ar.model(m).tUnits{1}];
            ar.model(m).vUnits{end,2} = [ar.model(m).xUnits{ix,2} '/' ar.model(m).tUnits{2}];
            ar.model(m).vUnits{end,3} = [ar.model(m).xUnits{ix,3} '/' ar.model(m).tUnits{3}];
            ar.model(m).v{end+1} = cell2mat(C{2});
            
            % setup N
            ar.model(m).N(1:length(ar.model(m).x),vcount) = 0;
            for jj=1:length(source)
                for j=find(ismember(ar.model(m).x, source{jj})) %R2013a compatible
                    ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) + sourceCoeffs(jj);
                end
            end
            for jj=1:length(target)
                for j=find(ismember(ar.model(m).x, target{jj})) %R2013a compatible
                    ar.model(m).N(j, vcount) = ar.model(m).N(j, vcount) - targetCoeffs(jj);
                end
            end
            
            % check for inconsistent educt compartments
            if(~isempty(ar.model(m).c))
                for j=1:size(ar.model(m).N,2)
                    if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)>0)))>1)
                        arParsingError( fid, 'efflux from different compartments in reaction %s', ...
                            ar.model(m).fv{end});
                    end
                    if(length(unique(ar.model(m).cLink(ar.model(m).N(:,j)<0)))>1)
                        arParsingError( fid, 'influx from different compartments in reaction %s', ...
                            ar.model(m).fv{end});
                    end
                end
            end
            
            vcount = vcount + 1;
        end
        
        %str = textscan(fid, '%s',1, 'CommentStyle', ar.config.comment_string);
        % No string left? Grab a new one
        if isempty( remainder )
            [ line, remainder, fid ] = readLine( fid, ar.config.comment_string );
        end
        [str, remainder] = grabtoken( remainder, '%s', 1 );
    end
elseif(strcmp(C{1},'ODES'))
    ar.model(m).isReactionBased = false;
    [str, fid] = arTextScan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
    ode_count = 0;
    while(~strcmp(str{1},'INVARIANTS') && ~strcmp(str{1},'DERIVED'))
        if(~strcmp(str{1},''))
            ode_count = ode_count + 1;
            ar.model(m).fv{end+1,1} = cell2mat(str{1});
            ar.model(m).fv_ma_reverse_pbasename{end+1} = '';
            ar.model(m).vUnits{end+1,1} = [ar.model(m).xUnits{ode_count,1} '/' ar.model(m).tUnits{1}];
            ar.model(m).vUnits{end,2} = [ar.model(m).xUnits{ode_count,2} '/' ar.model(m).tUnits{2}];
            ar.model(m).vUnits{end,3} = [ar.model(m).xUnits{ode_count,3} '/' ar.model(m).tUnits{3}];
        end
        [str, fid] = arTextScan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
    end
    if(ode_count ~= length(ar.model(m).x))
        arParsingError( fid, 'number of ODES ~= number of variables');
    end
    ar.model(m).N = eye(length(ar.model(m).x));
end
ar.model(m).qPlotV = ones(1,length(ar.model(m).fv));
if(isempty(ar.model(m).fv))
    ar.model(m).isReactionBased = false;
    ar.model(m).fv{end+1} = '0';
    ar.model(m).reversible(end+1) = 0;
    ar.model(m).fv_ma_reverse_pbasename{end+1} = '';
    ar.model(m).fv_source{end+1} = {};
    ar.model(m).fv_target{end+1} = {};
    ar.model(m).fv_sourceCoeffs{end+1} = [];
    ar.model(m).fv_targetCoeffs{end+1} = [];
    ar.model(m).vUnits{end+1,1} = '-';
    ar.model(m).vUnits{end,2}   = '-';
    ar.model(m).vUnits{end,3}   = '-';
    
    ar.model(m).N(1:length(ar.model(m).x),1) = 0;
end

% dynamic parameters
ar.model(m).pvs = cell(size(ar.model(m).fv));
ar.model(m).pv = {};
for jv=1:length(ar.model(m).fv)
    varlist = symvar(ar.model(m).fv{jv});
    ar.model(m).pvs{jv} = setdiff(varlist, union(ar.model(m).t, union(ar.model(m).x, ar.model(m).u))); %R2013a compatible
    ar.model(m).pv = union(ar.model(m).pv, ar.model(m).pvs{jv});
end
ar.model(m).px = union(union(ar.model(m).pv, ar.model(m).px), ar.model(m).px0); %R2013a compatible
ar.model(m).p = union(ar.model(m).px, ar.model(m).pu); %R2013a compatible

% setup rhs
C = cell(size(ar.model(m).N));
C_par = cell(size(ar.model(m).N));
if(length(ar.model(m).c)>1)    
    if(~isfield(ar.model(m),'isAmountBased') || ~ar.model(m).isAmountBased)
        for j=1:size(ar.model(m).N,1) % for every species j
            qinfluxwitheducts = ar.model(m).N(j,:) > 0 & sum(ar.model(m).N < 0,1) > 0;
            eductcompartment = zeros(size(qinfluxwitheducts));
            for jj=find(qinfluxwitheducts)
				eductcompartment(jj) = unique(ar.model(m).cLink(ar.model(m).N(:,jj)<0)); %R2013a compatible
            end
            
            cfaktor = cell(size(qinfluxwitheducts));
            cfaktor_par = cell(size(qinfluxwitheducts));
            for jj=1:size(ar.model(m).N,2) % for every reaction jj
                if(qinfluxwitheducts(jj) && eductcompartment(jj)~=ar.model(m).cLink(j))
                    cfaktor{jj} = [ar.model(m).pc{eductcompartment(jj)} '/' ...
                        ar.model(m).pc{ar.model(m).cLink(j)}];
                    cfaktor_par{jj} = [ar.model(m).c{eductcompartment(jj)} '/' ...
                        ar.model(m).c{ar.model(m).cLink(j)}];
                else
                    cfaktor{jj} = '1';
                    cfaktor_par{jj} = '1';
                end
            end
            C(j,:) = transpose(cfaktor);
            C_par(j,:) = transpose(cfaktor_par);
        end
    else
        for j=1:size(ar.model(m).N,1) % for every species j
            for jj=1:size(ar.model(m).N,2) % for every reaction jj
                C{j,jj} = ['1/' ar.model(m).pc{ar.model(m).cLink(j)}];
                C_par{j,jj} = ['1/' ar.model(m).c{ar.model(m).cLink(j)}];
            end
        end
    end
else
    for j=1:size(ar.model(m).N,1) % for every species j
        for jj=1:size(ar.model(m).N,2) % for every reaction jj
            C{j,jj} = '1';
            C_par{j,jj} = '1';
        end
    end
end
ar.model(m).fx = cell(length(ar.model(m).x),1);
ar.model(m).fx_par = cell(length(ar.model(m).x),1);

% initialize symbolic variables
% Joep: I have removed this. Not sure if this serves any purpose? Removing
% it seems to not affect any of the integration tests.
%if(~isempty(ar.model(m).x))
%    eval(['syms ' sprintf('%s ',ar.model(m).x{:})]);
%end
%if(~isempty(ar.model(m).p))
%    eval(['syms ' sprintf('%s ',ar.model(m).p{:})]);
%end
%if(~isempty(ar.model(m).u))
%    eval(['syms ' sprintf('%s ',ar.model(m).u{:})]);
%end
ar.model(m).Cm = C;
ar.model(m).Cm_par = C_par;
try
    fvSym = arMyStr2Sym(ar.model(m).fv);
catch
    error('Invalid expression in  expression equations:\n%s\n', sprintf('%s\n', ar.model(m).fv{:}));
end
    
tmpfx = (arMyStr2Sym(ar.model(m).N).* arMyStr2Sym(C)) * fvSym;
tmpfx_par = (arMyStr2Sym(ar.model(m).N).* arMyStr2Sym(C_par)) * fvSym;

for j=1:length(ar.model(m).x) % for every species j
    if ~isempty(tmpfx)
        ar.model(m).fx{j} = char(tmpfx(j));
        ar.model(m).fx_par{j} = char(tmpfx_par(j));
    else
        ar.model(m).fx{j} = char('0');
        ar.model(m).fx_par{j} = char('0');
    end
end


% DERIVED (previously INVARIANTS)
if(strcmp(str{1},'INVARIANTS'))
    arParsingError( fid, ['Section INVARIANTS in model definition file is deprecated! ' ...
        'Please replace by DERIVED and see usage in: ' ...
        'https://github.com/Data2Dynamics/d2d/wiki/Setting%20up%20models']);
end
derivedVariablesInRates = 0;
ar.model(m).z = {};
ar.model(m).zUnits = {};
ar.model(m).fz = {};
[C, fid] = arTextScan(fid, '%s %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'CONDITIONS') && ~strcmp(C{1},'SUBSTITUTIONS') && ~strcmp(C{1},'OBSERVABLES'))
    if(~strcmp(C{1},''))
        arValidateInput( C, 'derived', 'unique identifier', 'unit type (i.e. C)', 'unit (i.e. "units/cell")', 'plain text label for plots ("conc.")', 'derived function expression' );
        if(sum(ismember(ar.model(m).x, C{1}))>0) %R2013a compatible
            arParsingError( fid, 'derived variable %s already defined in STATES', cell2mat(C{1}));
        end
        if(sum(ismember(ar.model(m).u, C{1}))>0) %R2013a compatible
            arParsingError( fid, 'derived variable %s already defined in INPUTS', cell2mat(C{1}));
        end
        % Found derived variable in parameter list. See if it is used in
        % either the inputs or the reaction equations. If the former, fail
        % the loading, since we cannot resolve the state values at the time
        % the inputs are computed. If the latter => substitute them in!
        if(sum(ismember(ar.model(m).p, C{1}))>0) %R2013a compatible
            ar.model(m).p(ismember(ar.model(m).p, C{1})) = [];
            derivedVariablesInRates = 1;
            fail = 0;
            inputVariables = cellfun(@symvar, ar.model(m).fu, 'UniformOutput', false);
            for ju = 1 : length( inputVariables )
                fail = fail | max( ismember( inputVariables{ju}, C{1} ) );
            end
            if ( fail )
                arParsingError( fid, 'derived variable %s already defined as parameter in INPUT section', cell2mat(C{1}));
            end
        end
        ar.model(m).z(end+1) = C{1};
        ar.model(m).zUnits(end+1,1) = C{2};
        ar.model(m).zUnits(end,2) = C{3};
        ar.model(m).zUnits(end,3) = C{4};
        ar.model(m).fz(end+1,1) = C{5};
    end
    [C, fid] = arTextScan(fid, '%s %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
end

% Perform (repeated) self-substitutions
if numel( ar.model(m).fz ) > 0
    for a = 1 : numel( ar.model(m).fz )
        ar.model(m).fz{a} = char( arSubsRepeated(arMyStr2Sym(ar.model(m).fz{a}), ar.model(m).z, ar.model(m).fz, matVer.Version) );
    end
    arFprintf(2, '=> Self-substituting derived variables.\n' );
end

% Perform (repeated) derived substitutions
if ( derivedVariablesInRates )
    for a = 1 : length( ar.model(m).fv )
        ar.model(m).fv{a} = char( arSubsRepeated(arMyStr2Sym(ar.model(m).fv{a}), ar.model(m).z, ar.model(m).fz, matVer.Version) );
    end
    arFprintf(2, '=> Substituting derived variables in reaction equation.\n' );
end

% Simplify equations
%if ( simplifyEquations )
    
%end


ar.model(m).qPlotZ = ones(size(ar.model(m).z));

% derived variables parameters
varlist = cellfun(@symvar, ar.model(m).fz, 'UniformOutput', false);
ar.model(m).pz = setdiff(setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}), union(ar.model(m).x, union(ar.model(m).u, ar.model(m).z))); %R2013a compatible
ar.model(m).px = union(ar.model(m).px, ar.model(m).pz); %R2013a compatible
ar.model(m).p = union(ar.model(m).p, ar.model(m).pz); %R2013a compatible


% OBSERVABLES
% The following initialization should be done in any case. Otherwise
% problems occur if a model without obeservables is followed by a model
% with observables.
ar.model(m).y = {};
ar.model(m).yNames = {};
ar.model(m).yUnits = {};
ar.model(m).normalize = [];
ar.model(m).logfitting = [];
ar.model(m).logplotting = [];
ar.model(m).fy = {};
ar.model(m).fystd = {};

if(strcmp(C{1},'OBSERVABLES'))    
    [C, fid] = arTextScan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', ar.config.comment_string);
    while(~strcmp(C{1},'ERRORS'))
        if ( strcmp( C{1}, 'CONDITIONS' ) || strcmp( C{1}, 'SUBSTITUTIONS' ) )
            arParsingError( fid,  'When OBSERVABLES section is specified; ERRORS section must also be specified.' );
        end
        arValidateInput( C, 'observable', 'unique identifier', 'unit type (i.e. C)', 'unit (i.e. "units/cell")', 'plain text label for plots ("conc.")', 'indicator whether data should be scaled to 1 (0 or 1)', 'Indicator whether date should be treated in log-space (0 or 1)', 'Mathematical expression for observable' );
        ar.model(m).y(end+1) = C{1};
        ar.model(m).yUnits(end+1,1) = C{2};
        ar.model(m).yUnits(end,2) = C{3};
        ar.model(m).yUnits(end,3) = C{4};
        ar.model(m).normalize(end+1) = C{5};
        ar.model(m).logfitting(end+1) = C{6};
        ar.model(m).logplotting(end+1) = C{6};
        ar.model(m).fy(end+1,1) = C{7};
        if(~isempty(cell2mat(C{8})))
            ar.model(m).yNames(end+1) = C{8};
        else
            ar.model(m).yNames(end+1) = ar.model(m).y(end);
        end
        [C, fid] = arTextScan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', ar.config.comment_string);
        if(sum(ismember(ar.model(m).x, ar.model(m).y{end}))>0) %R2013a compatible
            arParsingError( fid, '%s already defined in STATES', ar.model(m).y{end});
        end
        if(sum(ismember(ar.model(m).u, ar.model(m).y{end}))>0) %R2013a compatible
            arParsingError( fid, '%s already defined in INPUTS', ar.model(m).y{end});
        end
        if(sum(ismember(ar.model(m).z, ar.model(m).y{end}))>0) %R2013a compatible
            arParsingError( fid, '%s already defined in DERIVED', ar.model(m).y{end});
        end
        if(sum(ismember(ar.model(m).p, ar.model(m).y{end}))>0) %R2013a compatible
            arParsingError( fid, '%s already defined as parameter', ar.model(m).y{end});
        end
    end
    
    % observation parameters
    varlist = cellfun(@symvar, ar.model(m).fy, 'UniformOutput', false);
    ar.model(m).py = setdiff(setdiff(vertcat(varlist{:}), union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
    
    % ERRORS
    ar.model(m).fystd = cell(size(ar.model(m).fy));
    [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
    while(~(strcmp(C{1},'CONDITIONS') || strcmp(C{1},'SUBSTITUTIONS')))
        qy = ismember(ar.model(m).y, C{1}); %R2013a compatible
        
        if(sum(qy)<1)
            arParsingError( fid, 'Unknown observable in error specification %s.', cell2mat(C{1}));
        elseif sum(qy)>1
            arParsingError( fid, 'Observable %s seems to occur more than once.', cell2mat(C{1}));            
        end
        
        %check and error if observable in log and fystd = rel + abs error
        y_var_name = setdiff(symvar(ar.model(m).fy{qy}),ar.model(m).py);
        reg_string = ['((?<=\W)|^)(',C{1}{1},'|'];
        for jreg = 1:length(y_var_name)
            if(jreg<length(y_var_name))
                reg_string = [reg_string ,y_var_name{jreg},'|'];
            else
                reg_string = [reg_string ,y_var_name{jreg},')'];
            end
        end
        reg_string = [reg_string '((?=\W)|$)'];
        if(~isempty(regexp(C{2}{1},reg_string,'ONCE')) && ar.model(m).logfitting(qy))
            warning(['You are trying to set up a relative error model within a log transformation. \n%s' ...
            'Comment out this error if you want to proceed anyway. To implement an absolute error in log, \n' ...
            'you can try the approach: \nyObs = sd_yObs + 1/2 * (a+sqrt((a)^2)), a = (offset - yObs-sd_yObs) \n, with hard set or fitted offset (on log-scale) \n'],C{2}{1})
            arParsingError( fid, 'Revise error model')
        end
        arValidateInput( C, 'error', 'observable identifier', 'expression for the error model' );
        
        ar.model(m).fystd(qy) = C{2};
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        if ( isempty( C{1} ) )
            arParsingError( fid,  'Missing field CONDITIONS' );
        end
    end
    
    if (length(ar.model(m).fystd)<length(ar.model(m).fy) || sum(cellfun(@isempty, ar.model(m).fystd))>0)
        diffErr = ar.model(m).y(cellfun(@isempty, ar.model(m).fystd)>0);
        if ( length(ar.model(m).fystd)<length(ar.model(m).fy) )
            diffErr = union( ar.model(m).y( length(ar.model(m).fystd) + 1 : end ), diffErr );
        end
        
        arParsingError( fid, 'Some observables do not have an error model defined. Observable(s) without error model: %s\n', sprintf( '%s ', diffErr{:} ) );
    end
    
    % error parameters
    varlist = cellfun(@symvar, ar.model(m).fystd, 'UniformOutput', false);
    ar.model(m).pystd = setdiff(vertcat(varlist{:}), union(union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ... %R2013a compatible
        ar.model(m).z), ar.model(m).y), ar.model(m).t));
    
    % add to parameters needed for model
    ar.model(m).p = union(union(ar.model(m).p, ar.model(m).py), ar.model(m).pystd);
end

% SUBSTITUTIONS (beta)
substitutions = 0;
if ( strcmp(C{1},'SUBSTITUTIONS') )
    if(str2double(matVer.Version)>=8.4)
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
    else
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16);
    end
    
    % Substitutions
    fromSubs = {};
    toSubs = {};
    ismodelpar = [];

    % Fetch desired substitutions
    while(~isempty(C{1}) && ~strcmp(C{1},'CONDITIONS'))
        fromSubs(end+1)     = C{1}; %#ok<AGROW>
        toSubs(end+1)       = C{2}; %#ok<AGROW>
        ismodelpar(end+1)   = sum(ismember(ar.model(m).p, C{1})); %#ok<AGROW>

        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
        
        if ( ~strcmp(C{1},'CONDITIONS') )
            arValidateInput( C, 'substitution', 'substitution identifier', 'expression for the substitution' );
        end
    end

    if ( sum(ismodelpar) > 0 )
        s = sprintf( '%s\n', fromSubs{ismodelpar>0} );
        arParsingError( fid,  'Cannot substitute model parameters. These following parameters belong under CONDITIONS:\n%s', s );
    end

    % Perform selfsubstitutions
    if ( ~isempty(fromSubs) )
        substitutions = 1;
        try
            toSubs = arSubsRepeated( toSubs, fromSubs, toSubs, str2double(matVer.Version) );
        catch
            arParsingError( fid,  'Error parsing substitution. Check the following expression for errors:\n%s', sprintf('%s\n', toSubs{:}) );
        end
    end
end


% CONDITIONS
if(str2double(matVer.Version)>=8.4)
    [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
else
    [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16);
end
ar.model(m).fp = transpose(ar.model(m).p);

% Remove init of empty model
if numel( ar.model(m).x ) == 1 && strcmp( ar.model(m).x, NOSTATE )
    ar.model(m).fp{ strcmp( ar.model(m).p, ['init_' NOSTATE ] ) } = '0';
end

if ( substitutions )
	% Conditions
    from        = {};
    to          = {};
    ismodelpar  = [];
    
    % Fetch desired conditions
    while(~isempty(C{1}) && ~(strcmp(C{1},'PARAMETERS') || strcmp(C{1}, 'RANDOM')))
        arValidateInput( C, 'condition', 'model parameter', 'new expression' );
        from(end+1)         = C{1}; %#ok<AGROW>
        to(end+1)           = C{2}; %#ok<AGROW>
        ismodelpar(end+1)   = sum(ismember(ar.model(m).p, C{1})); %#ok<AGROW>

        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
    end
    
    % Extra conditions specified at compile time
    if ( opts.conditions )
        for a = 1 : length( opts.conditions_args )
            from(end+1)         = opts.conditions_args{a,1}; %#ok<AGROW>
            to(end+1)           = opts.conditions_args{a,2}; %#ok<AGROW>
            ismodelpar(end+1)   = sum(ismember(ar.model(m).p, from(end))); %#ok<AGROW>
        end
    end
        
    % Perform substitutions (self-substitutions were already done)
    to = arSubsRepeated( to, fromSubs, toSubs, str2double(matVer.Version) );
    
    % Store substitutions in ar structure
    for a = 1 : length( from )
        qcondpara = ismember(ar.model(m).p, from{a}); %R2013a compatible
        if(sum(qcondpara)>0)
            ar.model(m).fp{qcondpara} = ['(' to{a} ')'];
        else
            warning('unknown parameter in conditions: %s (did you mean to place it under SUBSTITUTIONS?)', from{a}); %#ok<WNTAG>
        end
    end
else
    % Old code path
    while(~isempty(C{1}) && ~(strcmp(C{1},'PARAMETERS') || strcmp(C{1}, 'RANDOM')))
        arValidateInput( C, 'condition', 'model parameter', 'new expression' );
        qcondpara = ismember(ar.model(m).p, C{1}); %R2013a compatible
        if(sum(qcondpara)>0)
            ar.model(m).fp{qcondpara} = ['(' cell2mat(C{2}) ')'];
        else
            warning('unknown parameter in conditions: %s (did you mean to place it under SUBSTITUTIONS?)', cell2mat(C{1})); %#ok<WNTAG>
        end
        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
    end
end

ar.model(m).prand = {};
ar.model(m).rand_type = [];
if ( strcmp(C{1}, 'RANDOM' ) )    
    [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
    while(~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS'))
        ar.model(m).prand{end+1} = cell2mat(C{1});
        if(strcmp(C{2}, 'INDEPENDENT'))
            ar.model(m).rand_type(end+1) = 0;
        elseif(strcmp(C{2}, 'NORMAL'))
            ar.model(m).rand_type(end+1) = 1;
        else
            warning('unknown random type %s', cell2mat(C{2}));  %#ok<WNTAG>
        end
        [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
    end    
end
    
% extra conditional parameters
varlist = cellfun(@symvar, ar.model(m).fp, 'UniformOutput', false);
ar.model(m).pcond = setdiff(setdiff(setdiff(vertcat(varlist{:}), ar.model(m).p), ar.model(m).x), ar.model(m).u); %R2013a compatible

% PARAMETERS
if(~isfield(ar, 'pExternLabels'))
    ar.pExternLabels = {};
    ar.pExtern = [];
    ar.qFitExtern = [];
    ar.qLog10Extern = [];
    ar.lbExtern = [];
    ar.ubExtern = [];
end
[C, fid] = arTextScan(fid, '%s %f %n %n %f %f\n',1, 'CommentStyle', ar.config.comment_string);

while(~isempty(C{1}))
    ar.pExternLabels(end+1) = C{1};
    ar.pExtern(end+1) = C{2};
    ar.qFitExtern(end+1) = C{3};
    ar.qLog10Extern(end+1) = C{4};
    ar.lbExtern(end+1) = C{5};
    ar.ubExtern(end+1) = C{6};
    [C, fid] = arTextScan(fid, '%s %f %n %n %f %f\n',1, 'CommentStyle', ar.config.comment_string);
end

if ~isstruct( fid )
    fclose(fid);
end

% Check whether the user specified any variables with reserved words. This
% would be problematic later.
for a = 1 : length( ar.model(m).fu )
    arCheckReservedWords( symvar(ar.model(m).fu{a}), 'input function', ar.model(m).u{a} );
end
if ( isfield( ar.model(m), 'fy' ) )
    for a = 1 : length( ar.model(m).fy )
        arCheckReservedWords( symvar(ar.model(m).fy{a}), 'observation function', ar.model(m).y{a} );
    end
    for a = 1 : length( ar.model(m).fystd )
        arCheckReservedWords( symvar(ar.model(m).fystd{a}), 'observation standard deviation', ar.model(m).y{a} );
    end
end
if ( isfield( ar.model(m), 'fp' ) )
    for a = 1 : length( ar.model(m).fp )
        arCheckReservedWords( symvar(ar.model(m).fp{a}), 'parameter transformation', ar.model(m).p{a} );
    end
end
for a = 1 : length( ar.model(m).fx )
    arCheckReservedWords( symvar(ar.model(m).fx{a}), 'right hand side', ar.model(m).x{a} );
end

arCheckReservedWords( ar.model(m).p, 'parameters' );
arCheckReservedWords( ar.model(m).x, 'state variables' );
if ( isfield( ar.model(m), 'y' ) )
    arCheckReservedWords( ar.model(m).y, 'observables' );
end
if ( isfield( ar.model(m), 'z' ) )
    arCheckReservedWords( ar.model(m).z, 'derived variables' );
end
if ( isfield( ar.model(m), 'u' ) )
    arCheckReservedWords( ar.model(m).u, 'inputs' );
end
arCheckReservedWords( ar.model(m).c, 'compartments' );

ar = orderfields(ar);
ar.model = orderfields(ar.model);

if ( ar.config.checkForNegFluxes==2 )
    arCheckForNegFluxes(m, fid);
end

function [ line, remainder, fid ] = readLine( fid, commentStyle )
    line = '';
    while ( isempty(line) )
        [line, fid] = arTextScan(fid, '%[^\n]' ); 
        line = strtrim(line{1}{1});
        Q = strfind( line, commentStyle );
        if ( ~isempty(Q) )
            line = line(1:Q-1);
        end
    end
    remainder = line;
    
    
function [str, remainder] = grabtoken( inputString, varargin )
    if ( isempty( inputString ) )
        str{1} = {};
        remainder = '';
        return;
    end

	[str, pos] = textscan(inputString, varargin{:});
	remainder = inputString(pos+1:end);    

function num = checkNum( num, defaultValue )
    if ( ~isnumeric( num ) || isempty( num ) || isnan( num ) )
        num = defaultValue;
    end