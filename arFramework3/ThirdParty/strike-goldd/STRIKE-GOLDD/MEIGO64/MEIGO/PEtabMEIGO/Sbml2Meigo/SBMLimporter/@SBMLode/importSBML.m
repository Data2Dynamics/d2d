function importSBML(obj, model, verbose)
% importSBML parses information from the SBML definition and populates
% the SBMLode object from this information.
%
% Parameters:
%  filename: target name of the model
%
% Return values:
%  void

if nargin == 2
    verbose = false;
end

if(isfield(model,'fbc_objective'))
    error('Flux Balance Constraints are currently not supported')
end

%% COMPARTMENTS

if verbose
    fprintf('loading compartments ...\n')
end

% initialize
compartments_sym = sym({model.compartment.id});
obj.compartment = compartments_sym;

% set initial assignments
if(isfield(model,'initialAssignment'))
    initassignments_sym = sym({model.initialAssignment.symbol});
    initassignments_math = cleanedsym({model.initialAssignment.math});
else
    initassignments_sym = sym([]);
    initassignments_math = sym([]);
end
setInitialAssignment(obj,model,'compartment',initassignments_sym,initassignments_math)

% % extract compartment sizes, default to 1
% obj.compartment = subs(obj.compartment,compartments_sym(logical([model.compartment.isSetSize])),sym([model.compartment.size]));
%
% % set remaining ones to default value 1
% obj.compartment = subs(obj.compartment,compartments_sym,sym(ones(size(compartments_sym))));

%% RULES
if verbose
    fprintf('applying rules ...\n')
end

rule_types = {model.rule.typecode};
if(any(strcmp(rule_types,'SBML_ALGEBRAIC_RULE')))
    %DAE part TBD
    error('Algebraic rules are currently not supported!');
end

all_rulevars = sym({model.rule.variable});
if(any(arrayfun(@(x) ismember(x,compartments_sym),all_rulevars)))
    error('Rules for compartments are currently not supported!');
end
all_rulemath = cleanedsym({model.rule.formula});
% remove rate rules
rulevars = all_rulevars(not(strcmp({model.rule.typecode},'SBML_RATE_RULE')));
rulemath = all_rulemath(not(strcmp({model.rule.typecode},'SBML_RATE_RULE')));
repeat_idx = ismember(rulevars,symvar(rulemath));
while(any(repeat_idx))
    rulemath= subs(rulemath,rulevars,rulemath);
    repeat_idx = ismember(rulevars,symvar(rulemath));
end

applyRule(obj,model,'compartment',rulevars,rulemath)

%% SPECIES
if verbose
    fprintf('loading species ...\n')
end

% extract species
species_sym = sym({model.species.id});
species_sym = species_sym(:);
obj.state = species_sym;

nx = length(obj.state);

% extract corresponding volumes
compartments = sym(sanitizeString({model.species.compartment}));
obj.volume = subs(compartments(:),sym({model.compartment.id}),obj.compartment);
if(any(arrayfun(@(x) ~isempty(regexp(char(x),'[\w]+\(')),obj.volume)))
    error('Functions in volume definitions is currently not supported.')
end

initConcentration = [model.species.initialConcentration];
initConcentration = initConcentration(:);
initAmount = [model.species.initialAmount];
initAmount = initAmount(:);
obj.initState = species_sym;
hasAssignmentRule = ismember(obj.state,all_rulevars(strcmp({model.rule.typecode},'SBML_ASSIGNMENT_RULE')))';
hasRateRule = ismember(obj.state,all_rulevars(strcmp({model.rule.typecode},'SBML_RATE_RULE')))';

% set initial assignments
setInitialAssignment(obj,model,'initState',initassignments_sym,initassignments_math)



% remove conditions species (boundary condition + no initialisation)
if(~isempty(obj.state))
    cond_idx = all([logical([model.species.boundaryCondition]);transpose(logical(obj.state==obj.initState));~hasAssignmentRule;~hasRateRule]);
    bound_idx = all([logical([model.species.boundaryCondition]);transpose(logical(obj.state~=obj.initState));~hasAssignmentRule;~hasRateRule]);
    condition_sym = obj.state(cond_idx);
    conditions = obj.state(cond_idx);
    boundary_sym = obj.state(bound_idx);
    boundaries = obj.initState(bound_idx);
else
    
    cond_idx = [];
    bound_idx = [];
    condition_sym = sym([]);
    conditions = sym([]);
    boundary_sym = sym([]);
    boundaries = sym([]);
end

nk = length(conditions);

applyRule(obj,model,'condition',rulevars,rulemath)

% set initial concentrations
concentration_idx = logical([model.species.isSetInitialConcentration]);
onlysubstance_idx = logical([model.species.hasOnlySubstanceUnits]);
obj.initState = subs(obj.initState,species_sym(concentration_idx),sym(initConcentration(concentration_idx)));

% set initial amounts
amount_idx = logical([model.species.isSetInitialAmount]);
obj.initState = subs(obj.initState,species_sym(amount_idx),sym(initAmount(amount_idx))./obj.volume(amount_idx));
% this.initState = subs(this.initState,species_sym(amount_idx),sym(initAmount(amount_idx)));

% apply rules
applyRule(obj,model,'initState',rulevars,rulemath)

while(any(ismember(symvar(obj.initState),obj.state)))
    obj.initState = subs(obj.initState,obj.state,obj.initState);
end
while(any(ismember(symvar(boundaries),obj.state)))
    boundaries = subs(boundaries,obj.state,obj.initState);
end

%% PARAMETERS
if verbose
    fprintf('loading parameters ...\n')
end

% extract this.param
parameter_sym = sym({model.parameter.id});
parameter_val = transpose([model.parameter.value]);
parameter_sym = parameter_sym(:);
obj.param = parameter_sym;

np = length(obj.param);

%% CONSTANTS

% remove constant species
const_idx = logical([model.species.constant]) & not(cond_idx);
constant_sym = obj.state(const_idx);

obj.knom = double(subs([obj.initState(cond_idx);obj.initState(const_idx)],parameter_sym,parameter_val));





%% REACTIONS
if verbose
    fprintf('parsing reactions ...\n')
end

nr = length(model.reaction);
if(nr>0)
    for ir = 1:nr
        if(model.reaction(ir).isSetFast)
            if(model.reaction(ir).fast)
                error('Fast reactions are currently not supported!');
            end
        end
    end
end

kLaw = [cellfun(@(x) x.math,{model.reaction.kineticLaw},'UniformOutput',false)];

checkIllegalFunctions(kLaw);

obj.flux = cleanedsym(kLaw);
obj.flux = obj.flux(:);
% add local parameters to global parameters, make them global by
% extending them by the reaction_id string
species_idx = transpose(sym(1:nx));
if(length({model.reaction.id})>0)
    try
        tmp = cellfun(@(x,y) sym(cellfun(@(x) [x '_' y], ...
            {x.parameter.id}, ...
            'UniformOutput',false)), ...
            {model.reaction.kineticLaw}, ...
            {model.reaction.id}, ...
            'UniformOutput',false);
        plocal = transpose([tmp{:}]);
        tmp = cellfun(@(x) cellfun(@double,{x.parameter.value}),{model.reaction.kineticLaw},'UniformOutput',false);
        pvallocal = transpose([tmp{:}]);
        % replace local parameters by globalized ones
        tmp = cellfun(@(x,y,z) subs(x,sym({y.parameter.id}), ...
            sym(cellfun(@(x) [x '_' z],{y.parameter.id},'UniformOutput',false))),...
            transpose(num2cell(obj.flux)),...
            {model.reaction.kineticLaw},...
            {model.reaction.id},...
            'UniformOutput',false);
        obj.flux = [tmp{:}];
        obj.flux = obj.flux(:);
        
    catch
        tmp = cellfun(@(x,y) sym(cellfun(@(x) [x '_' y],{x.localParameter.id},'UniformOutput',false)),{model.reaction.kineticLaw},arrayfun(@(x) ['r' num2str(x)],1:length({model.reaction.id}),'UniformOutput',false),'UniformOutput',false);
        plocal = transpose([tmp{:}]);
        tmp = cellfun(@(x) cellfun(@double,{x.localParameter.value}),{model.reaction.kineticLaw},'UniformOutput',false);
        pvallocal = transpose([tmp{:}]);
        % replace local parameters by globalized ones
        tmp = cellfun(@(x,y,z) subs(x,sym({y.localParameter.id}),sym(cellfun(@(x) [x '_' z],{y.localParameter.id},'UniformOutput',false))),transpose(num2cell(obj.flux)),{model.reaction.kineticLaw},arrayfun(@(x) ['r' num2str(x)],1:length({model.reaction.id}),'UniformOutput',false),'UniformOutput',false);
        obj.flux = [tmp{:}];
        obj.flux = obj.flux(:);
        
    end
    
    obj.param = [obj.param;plocal];
    parameter_sym = [parameter_sym;plocal];
    parameter_val = [parameter_val;pvallocal];
    np = length(obj.param);
    
    reactants = cellfun(@(x) {x.species},{model.reaction.reactant},'UniformOutput',false);
    % species index of the reactant
    reactant_sidx = double(subs(sym(cat(2,reactants{:})),species_sym,species_idx));
    % reaction index
    tmp = cumsum(cell2mat(cellfun(@(x) [ones(1,1),zeros(1,max(length(x)-1,0))],reactants,'UniformOutput',false)));
    wreact = cell2mat(cellfun(@(x) [ones(1,length(x)),zeros(1,isempty(x))],reactants,'UniformOutput',false));
    reactant_ridx = tmp(logical(wreact));
    products = cellfun(@(x) {x.species},{model.reaction.product},'UniformOutput',false);
    % species index of the product
    product_sidx = double(subs(sym(cat(2,products{:})),species_sym,species_idx));
    % reaction index
    tmp = cumsum(cell2mat(cellfun(@(x) [ones(1,1),zeros(1,max(length(x)-1,0))],products,'UniformOutput',false)));
    wprod = cell2mat(cellfun(@(x) [ones(1,length(x)),zeros(1,isempty(x))],products,'UniformOutput',false));
    product_ridx = tmp(logical(wprod));
    if(model.SBML_level>=3)
        reactant_stochiometry = cellfun(@(x) stoich_initAssign_rule(x,initassignments_sym,initassignments_math,rulevars,rulemath),{model.reaction.reactant},'UniformOutput',false);
        %         reactant_math = cellfun(@(x) sym({x.stoichiometry}),{model.reaction.reactant},'UniformOutput',false);
        reactant_id = cellfun(@getId,{model.reaction.reactant},'UniformOutput',false);
        product_stochiometry = cellfun(@(x) stoich_initAssign_rule(x,initassignments_sym,initassignments_math,rulevars,rulemath),{model.reaction.product},'UniformOutput',false);
        %         product_math = cellfun(@(x) sym({x.stoichiometry}),{model.reaction.product},'UniformOutput',false);
        product_id = cellfun(@getId,{model.reaction.product},'UniformOutput',false);
    else
        % addition is necessary due to 1x0 struct that is returned by libSBML which is not properly handled by MATLAB,
        % the concatenation is necessary because MATLAB treats 1x0 structs as empty input
        symbolic_expr = @(x) num2cell(cell2sym(cellfun(@(z) math_expr(z),arrayfun(@(y) y.stoichiometryMath,x,'UniformOutput',false),'UniformOutput',false)) + sym(arrayfun(@(y) y.stoichiometry,x)).*arrayfun(@(y) isempty(y.stoichiometryMath),x));
        reactant_stochiometry = cellfun(@(x) {symbolic_expr(x)},{model.reaction.reactant},'UniformOutput',false);
        reactant_id = cellfun(@getId,{model.reaction.reactant},'UniformOutput',false);
        product_stochiometry = cellfun(@(x) {symbolic_expr(x)},{model.reaction.product},'UniformOutput',false);
        product_id = cellfun(@getId,{model.reaction.product},'UniformOutput',false);
    end
    eS = sym(zeros(nx,nr));
    pS = sym(zeros(nx,nr));
    tmp_rs = cellfun(@(x)[x{:}],reactant_stochiometry,'UniformOutput',false);
    tmp_rs = [tmp_rs{:}];
    tmp_rid = cat(2,reactant_id{:});
    for iidx = 1:length(reactant_sidx)
        if(strcmp(tmp_rid{iidx},''))
            tmp = tmp_rs(iidx);
        else
            tmp = sym(tmp_rid{iidx});
        end
        eS(reactant_sidx(iidx),reactant_ridx(iidx)) = eS(reactant_sidx(iidx),reactant_ridx(iidx)) + tmp;
    end
    tmp_ps = cellfun(@(x)[x{:}],product_stochiometry,'UniformOutput',false);
    tmp_ps = [tmp_ps{:}];
    tmp_pid = cat(2,product_id{:});
    for iidx = 1:length(product_sidx)
        if(strcmp(tmp_pid{iidx},''))
            tmp = tmp_ps(iidx);
        else
            tmp = sym(tmp_pid{iidx});
        end
        pS(product_sidx(iidx),product_ridx(iidx)) = pS(product_sidx(iidx),product_ridx(iidx)) + tmp;
    end
    
    obj.stochiometry = - eS + pS;
    
    obj.xdot = obj.stochiometry*obj.flux;
else
    obj.xdot = sym(zeros(size(obj.state)));
end

reactionsymbols = sym({model.reaction.id}');

if(length({model.reaction.id})>0)
    stoichsymbols = [reactant_id{:},product_id{:}];
    stoichmath = [tmp_rs,tmp_ps];
    
    stoichidx = not(strcmp(stoichsymbols,''));
    stoichsymbols = stoichsymbols(stoichidx);
    stoichmath = stoichmath(stoichidx);
else
    stoichsymbols = sym([]);
    stoichmath = sym([]);
end

%% RATE RULES
if verbose
    fprintf('converting to concentrations ...\n')
end

%extract model conversion factor
if(isfield(model,'conversionFactor'))
    if(strcmp(model.conversionFactor,''))
        conversionfactor = sym(ones(nx,1));
    else
        conversionfactor = sym(model.conversionFactor)*ones(nx,1);
    end
else
    conversionfactor = ones(nx,1);
end

if(isfield(model.species,'conversionFactor'))
    if(any(not(strcmp({model.species.conversionFactor},''))))
        tmp = {model.species.conversionFactor};
        idx = ~strcmp({model.species.conversionFactor},'');
        conversionfactor(idx) = sym(tmp(idx));
    end
end

for irule = 1:length(model.rule)
    if(strcmp(model.rule(irule).typecode,'SBML_RATE_RULE'))
        state_rate_idx = find(obj.state == sym(model.rule(irule).variable));
        param_rate_idx = find(parameter_sym == sym(model.rule(irule).variable));
        stoich_rate_idx = find(stoichsymbols == sym(model.rule(irule).variable));
        if(~isempty(state_rate_idx))
            obj.xdot(state_rate_idx) = cleanedsym(model.rule(irule).formula);
            if(~onlysubstance_idx(state_rate_idx))
                obj.xdot(state_rate_idx) = obj.xdot(state_rate_idx).*obj.volume(state_rate_idx);
            end
        elseif(~isempty(param_rate_idx))
            obj.state = [obj.state; parameter_sym(param_rate_idx)];
            obj.xdot = [obj.xdot; cleanedsym(model.rule(irule).formula)];
            if(ismember(parameter_sym(param_rate_idx),initassignments_sym))
                obj.initState = [obj.initState; initassignments_math(find(initassignments_sym==parameter_sym(param_rate_idx)))];
            else
                obj.initState = [obj.initState; parameter_val(param_rate_idx)];
            end
            obj.volume = [obj.volume; 1];
            concentration_idx = [concentration_idx, false];
            onlysubstance_idx = [onlysubstance_idx, false];
            conversionfactor = [conversionfactor; 1];
            cond_idx = [cond_idx, false ];
            const_idx = [const_idx, false ];
            bound_idx = [bound_idx, false];
            nx = nx + 1;
            parameter_val(param_rate_idx) = [];
            parameter_sym(param_rate_idx) = [];
            obj.param(param_rate_idx) = [];
            np = np - 1;
            setInitialAssignment(obj,model,'initState',initassignments_sym,initassignments_math);
        elseif(~isempty(stoich_rate_idx))
            obj.state = [obj.state; stoichsymbols(stoich_rate_idx)];
            obj.xdot = [obj.xdot; cleanedsym(model.rule(irule).formula)];
            obj.initState = [obj.initState; stoichmath(stoich_rate_idx)];
            obj.volume = [obj.volume; 1];
            concentration_idx = [concentration_idx, false];
            onlysubstance_idx = [onlysubstance_idx, false];
            conversionfactor = [conversionfactor; 1];
            cond_idx = [cond_idx, false ];
            const_idx = [const_idx, false ];
            bound_idx = [bound_idx, false];
            nx = nx + 1;
            stoichmath(stoich_rate_idx) = [];
            stoichsymbols(stoich_rate_idx) = [];
        end
    end
    if(strcmp(model.rule(irule).typecode,'SBML_ASSIGNMENT_RULE'))
        state_rate_idx = find(obj.state == sym(model.rule(irule).variable));
        param_rate_idx = find(parameter_sym == sym(model.rule(irule).variable));
        if(~isempty(state_rate_idx))
            obj.state(state_rate_idx) = [];
            obj.xdot(state_rate_idx) = [];
            obj.initState(state_rate_idx) = [];
            obj.volume(state_rate_idx) = [];
            concentration_idx(state_rate_idx) = [];
            onlysubstance_idx(state_rate_idx) = [];
            conversionfactor(state_rate_idx) = [];
            cond_idx(state_rate_idx) = [];
            const_idx(state_rate_idx) = [];
            bound_idx(state_rate_idx) = [];
            nx = nx-1;
            obj.observable = [obj.observable;cleanedsym(model.rule(irule).formula)];
            obj.observable_name = [obj.observable_name;sym(model.rule(irule).variable)];
        end
    end
end

applyRule(obj,model,'xdot',rulevars,rulemath)


%% CONVERSION FACTORS/VOLUMES

% this.xdot = conversionfactor.*subs(this.xdot,this.state,this.state.*this.volume)./this.volume;
obj.xdot = conversionfactor.*subs(obj.xdot,obj.state(onlysubstance_idx),obj.state(onlysubstance_idx).*obj.volume(onlysubstance_idx))./obj.volume;

% this.xdot = conversionfactor.*this.xdot;

%% EVENTS
if verbose
    fprintf('loading events ...\n')
end

if(sum(cellfun(@(x)numel(x),{model.event.delay})>0))
    error('Events with delays are currently not supported!');
end
if(~strcmp(model.delay_symbol,''))
    error('Delay symbols are currently not supported!');
end
if(model.SBML_level>=3)
    if(sum(cellfun(@(x)numel(x),{model.event.priority})>0))
        error('Event priorities are currently not supported!');
    end
end

try
    tmp = cellfun(@(x) sym(sanitizeString(x)),{model.event.trigger},'UniformOutput',false);
    obj.trigger = [tmp{:}];
catch
    tmp = cellfun(@(x) sym(sanitizeString(x.math)),{model.event.trigger},'UniformOutput',false);
    obj.trigger = [tmp{:}];
end
obj.trigger = obj.trigger(:);
obj.trigger = subs(obj.trigger,sym('ge'),sym('am_ge'));
obj.trigger = subs(obj.trigger,sym('gt'),sym('am_gt'));
obj.trigger = subs(obj.trigger,sym('le'),sym('am_le'));
obj.trigger = subs(obj.trigger,sym('lt'),sym('am_lt'));

obj.bolus = sym(zeros([length(obj.state),length(obj.trigger)]));
if(length(obj.trigger)>0)
    for ievent = 1:length(obj.trigger)
        tmp = cellfun(@(x) {x.variable},{model.event(ievent).eventAssignment},'UniformOutput',false);
        assignments = sym(cat(2,tmp{:}));
        
        tmp = cellfun(@(x) {x.math},{model.event(ievent).eventAssignment},'UniformOutput',false);
        assignments_math = cleanedsym(cat(2,tmp{:}));
        
        for iassign = 1:length(assignments)
            state_assign_idx = find(assignments(iassign)==obj.state);
            param_assign_idx = find(assignments(iassign)==obj.param);
            cond_assign_idx = find(assignments(iassign)==condition_sym);
            bound_assign_idx = find(assignments(iassign)==boundary_sym);
            stoich_assign_idx = find(assignments(iassign)==stoichsymbols);
            vol_assign_idx = find(assignments(iassign)==compartments_sym);
            
            if(np>0 && ~isempty(param_assign_idx))
                error('Assignments of parameters via events are currently not supported')
                obj.param(param_assign_idx) = obj.param(param_assign_idx)*heaviside(-obj.trigger(ievent)) + assignments_math(iassign)*heaviside(obj.trigger(ievent));
            end
            
            if(nk>0 && ~isempty(cond_assign_idx))
                error('Assignments of constants via events are currently not supported')
                conditions(cond_assign_idx) = conditions(cond_assign_idx)*heaviside(-obj.trigger(ievent)) + assignments_math(iassign)*heaviside(obj.trigger(ievent));
            end
            
            if(length(boundaries)>0 && ~isempty(bound_assign_idx))
                error('Assignments of boundary conditions via events are currently not supported')
                boundaries(bound_assign_idx) = conditions(bound_assign_idx)*heaviside(-obj.trigger(ievent)) + assignments_math(iassign)*heaviside(obj.trigger(ievent));
            end
            
            if(length(stoichsymbols)>0 && ~isempty(stoich_assign_idx))
                error('Assignments of stoichiometries via events are currently not supported')
                stoichmath(stoich_assign_idx) = stoichmath(stoich_assign_idx)*heaviside(-obj.trigger(ievent)) + assignments_math(iassign)*heaviside(obj.trigger(ievent));
            end
            
            if(length(compartments_sym)>0 && ~isempty(vol_assign_idx))
                error('Assignments of compartment volumes via events are currently not supported')
            end
            
            if(length(obj.state)>0 && ~isempty(state_assign_idx))
                
                obj.bolus(state_assign_idx,ievent) = -obj.state(state_assign_idx);
                addToBolus = sym(zeros(size(obj.bolus(:,ievent))));
                addToBolus(state_assign_idx) = assignments_math(iassign);
                
                obj.bolus(:,ievent) = obj.bolus(:,ievent) + addToBolus;
            end
            
        end
    end
else
    addToBolus = sym([]);
end




%% FUNCTIONS
if verbose
    fprintf('loading functions ...\n')
end

tmp = cellfun(@(x) x(8:end-1),{model.functionDefinition.math},'UniformOutput',false);
lambdas = cellfun(@(x)argScan(x),tmp);

if(~isempty(lambdas))
    obj.funmath = cellfun(@(x) x{end},lambdas,'UniformOutput',false);
    % temp replacement for any user defined function
    tmpfun = cellfun(@(x) ['fun_' num2str(x)],num2cell(1:length(model.functionDefinition)),'UniformOutput',false);
    obj.funmath = strrep(obj.funmath,{model.functionDefinition.id},tmpfun);
    % replace helper functions
    
    checkIllegalFunctions(obj.funmath);
    obj.funmath = replaceLogicalFunctions(obj.funmath);
    
    obj.funmath = strrep(obj.funmath,tmpfun,{model.functionDefinition.id});
    obj.funarg = cellfun(@(x,y) [y '(' strjoin(transpose(x(1:end-1)),',') ')'],lambdas,replaceReservedFunctionIDs({model.functionDefinition.id}),'UniformOutput',false);
    
    % make functions available in this file
    
    for ifun = 1:length(obj.funmath)
        token = regexp(obj.funarg(ifun),'\(([0-9\w\,]*)\)','tokens');
        start = regexp(obj.funarg(ifun),'\(([0-9\w\,]*)\)');
        eval([replaceReservedFunctions(obj.funarg{ifun}(1:(start{1}-1))) ' = @(' token{1}{1}{1} ')' obj.funmath{ifun} ';']);
    end
end


%% CLEAN-UP
if verbose
    fprintf('cleaning up ...\n')
end

% remove constant/condition states
obj.state(any([cond_idx;const_idx;bound_idx])) = [];
obj.kvolume = obj.volume(any([cond_idx;const_idx;bound_idx]));
obj.volume(any([cond_idx;const_idx;bound_idx])) = [];
obj.initState(any([cond_idx;const_idx;bound_idx])) = [];
obj.xdot(any([cond_idx;const_idx;bound_idx])) = [];
obj.bolus(any([cond_idx;const_idx;bound_idx]),:) = [];

% substitute with actual expressions, do this twice to resolve co-dependencies, do we need a loop here?
makeSubs(obj,boundary_sym,boundaries);
makeSubs(obj,condition_sym,conditions);
makeSubs(obj,compartments_sym,obj.compartment);
%makeSubs(this,stoichsymbols,stoichmath);
makeSubs(obj,reactionsymbols,obj.flux);

% set initial assignments
for iIA = 1:length(initassignments_sym)
    if(ismember(initassignments_sym(iIA),obj.param))
        if(ismember(sym(model.time_symbol),symvar(initassignments_math(iIA))))
            error('Time dependent initial assignments are currently not supported!')
        end
        param_idx =  find(initassignments_sym(iIA)==obj.param);
        parameter_sym(param_idx) = [];
        parameter_val(param_idx) = [];
        obj.param(param_idx) = [];
        obj.xdot = subs(obj.xdot,initassignments_sym(iIA),initassignments_math(iIA));
        obj.trigger = subs(obj.trigger,initassignments_sym(iIA),initassignments_math(iIA));
        obj.bolus = subs(obj.bolus,initassignments_sym(iIA),initassignments_math(iIA));
        obj.initState = subs(obj.initState,initassignments_sym(iIA),initassignments_math(iIA));
        obj.param = subs(obj.param,initassignments_sym(iIA),initassignments_math(iIA));
        rulemath = subs(rulemath,initassignments_sym(iIA),initassignments_math(iIA));
        np = np-1;
    end
end
applyRule(obj,model,'param',rulevars,rulemath)

makeSubs(obj,parameter_sym(1:np),obj.param);

% apply rules to dynamics
for irule = 1:length(rulevars)
    eval(['syms ' strjoin(arrayfun(@char,rulevars(irule),'UniformOutput',false),' ') ' ' strrep(strrep(strjoin(arrayfun(@char,symvar(sym(rulemath(irule))),'UniformOutput',false),' '),'true',''),'false','')]);
    eval(['rule = ' char(rulemath(irule)) ';']);
    rule_idx = find(rulevars(irule)==obj.state);
    if(nx>0)
        obj.xdot(rule_idx) = jacobian(rule,obj.state)*obj.xdot;
        obj.initState(rule_idx) = rulemath(irule);
        if(~isempty(obj.bolus))
            obj.bolus(rule_idx,:) = jacobian(rule,obj.state)*obj.bolus;
        end
    end
end

state_vars = [symvar(obj.xdot),symvar(obj.initState)];
event_vars = [symvar(obj.bolus),symvar(obj.trigger)];

applyRule(obj,model,'xdot',rulevars,rulemath)

isUsedParam = or(ismember(parameter_sym,event_vars),ismember(parameter_sym,state_vars));
isPartOfRule = and(ismember(parameter_sym,symvar(cleanedsym({model.rule.formula}))),ismember(parameter_sym,symvar(sym({model.rule.variable}))));
isRuleVar = ismember(parameter_sym,sym({model.rule.variable}));
hasAssignment = ismember(parameter_sym,initassignments_sym);
obj.parameter = parameter_sym(and(not(isRuleVar),not(isPartOfRule)));
obj.pnom = parameter_val(and(not(isRuleVar),not(isPartOfRule)));

obj.condition = [condition_sym;constant_sym];
obs_idx = all([isRuleVar,not(isPartOfRule),not(isUsedParam),not(hasAssignment)],2);
obj.observable = [obj.observable;obj.param(obs_idx(1:length(obj.param)))];
obj.observable_name = [obj.observable_name;parameter_sym(obs_idx(1:length(obj.param)))];

equal_idx = logical(obj.observable(:)==obj.observable_name(:)); % this is the unused stuff
obj.observable(equal_idx) = [];
obj.observable_name(equal_idx) = [];

obj.observable = subs(obj.observable,parameter_sym(1:np),obj.param);
obj.observable = subs(obj.observable,condition_sym,conditions);
obj.observable = subs(obj.observable,boundary_sym,boundaries);
obj.observable = subs(obj.observable,compartments_sym,obj.compartment);
obj.observable = subs(obj.observable,stoichsymbols,stoichmath);
obj.observable = subs(obj.observable,reactionsymbols,obj.flux);

applyRule(obj,model,'observable',rulevars,rulemath);

obj.time_symbol = model.time_symbol;

prohibited = sym({'null','beta'});
alt_prohib = sym({'null_sym','beta_sym'});
obj.parameter = subs(obj.parameter,prohibited,alt_prohib);
obj.state = subs(obj.state,prohibited,alt_prohib);
obj.condition = subs(obj.condition,prohibited,alt_prohib);
while(any(ismember(symvar(obj.initState),obj.state)))
    obj.initState = subs(obj.initState,obj.state,obj.initState);
end
end

function setInitialAssignment(this,~,field,initassignments_sym,initassignments_math)
this.(field) = subs(this.(field),initassignments_sym,initassignments_math);
end


function applyRule(this,~,field,rulevars,rulemath)
this.(field) = subs(this.(field),rulevars,rulemath);
end

function args = argScan(str)
brl = computeBracketLevel(str);

l1_idx = find(brl==0); % store indexes

% find commas but only in lowest bracket level
c_idx = strfind(str(brl==0),',');
str(l1_idx(c_idx)) = '@'; % replace by something that should never occur in equations
args = textscan(str,'%s','Whitespace','@');
end

function makeSubs(this,old,new)
this.xdot = subs(this.xdot,old,new);
this.trigger = subs(this.trigger,old,new);
this.bolus = subs(this.bolus,old,new);
this.initState = subs(this.initState,old,new);
end

function checkIllegalFunctions(str)

if(any(cell2mat(strfind(str,'factorial'))))
    error('Factorial functions are currently not supported!')
end
if(any(cell2mat(strfind(str,'ceil'))))
    error('Ceil functions are currently not supported!')
end
if(any(cell2mat(strfind(str,'floor'))))
    error('Floor functions are currently not supported!')
end
end

function csym = cleanedsym(str)
if(nargin>0)
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=9.4)
        csym = str2sym(sanitizeString(str));
    else
        csym = sym(sanitizeString(str));
    end
else
    csym = sym(0);
end
end

function str = sanitizeString(str)
% wrapper for replaceDiscontinuousFunctions() and replaceReservedFunctions()
str = replaceLogicalFunctions(str);
str = replaceReservedFunctions(str);
end

function str = replaceLogicalFunctions(str)
% replace imcompatible piecewise defintion
% execute twice for directly nested calls (overlapping regexp expressions)
for logicalf = {'piecewise','and','or','lt','gt','ge','le','ge','le','xor','eq'}
    str = regexprep(str,['^' logicalf{1} '('],['am_' logicalf{1} '(']);
    str = regexprep(str,['([\W]+)' logicalf{1} '('],['$1am_' logicalf{1} '(']);
    str = regexprep(str,['([\W]+)' logicalf{1} '('],['$1am_' logicalf{1} '(']);
end
end

function str = replaceReservedFunctions(str)
% replace reserved matlab functions

if(strcmp(str,'this'))
    error('SBML functions may not be called ''this''');
end

for logicalf = {'divide','minus','multiply','plus'}
    str = regexprep(str,['^' logicalf{1} '('],['am_' logicalf{1} '(']);
    % execute twice for directly nested calls (overlapping regexp expressions)
    str = regexprep(str,['([\W]+)' logicalf{1} '('],['$1am_' logicalf{1} '(']);
    str = regexprep(str,['([\W]+)' logicalf{1} '('],['$1am_' logicalf{1} '(']);
end
end

function str = replaceReservedFunctionIDs(str)
% replace reserved matlab functions
for logicalf = {'divide','minus','multiply','plus'}
    str = regexprep(str,['^' logicalf{1} '$'],['am_' logicalf{1} '']);
end
end

function z = delay(x,y)
error('Events with delays are currently not supported!');
end

function x = stoich_initAssign_rule(y,initassignments_sym,initassignments_math,rulevars,rulemath)
x = {};
for iy = 1:length(y)
    x{iy} = sym(y(iy).stoichiometry);
    if(~isempty(sym(y(iy).id)))
        if(~isempty(initassignments_sym))
            if(ismember(sym(y(iy).id),initassignments_sym))
                x{iy} = subs(sym(y(iy).id),initassignments_sym,initassignments_math);
            end
        end
        if(~isempty(rulevars))
            if(ismember(sym(y(iy).id),rulevars))
                x{iy} = subs(sym(y(iy).id),rulevars,rulemath);
            end
        end
    end
end
end

function expr = math_expr(y)
if(isfield(y,'math'))
    expr = cleanedsym(y.math);
else
    expr = cleanedsym();
end
end

function id = getId(x)
if(isfield(x,'id'))
    id = {x.id};
else
    id = {x.species};
end
end

