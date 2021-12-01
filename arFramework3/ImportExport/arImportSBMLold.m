% varargout = arImportSBMLold(filename, varargin)
%
% Import SBML model and translate to .def files
%
% Options which can be specified are:
%   tEnd                - Final simulation time (default = 100)
%   compartmentbyname   - Use compartment names as specified in the SBML rather than unique identifier (default = false)
%   overwrite           - Overwrite def file if it exists (default = false)
%   keepcompartments    - Reference compartments in reactions rather than numerical values (default = false)
%
% Optional Outputs:
%   ms   - struct with typical sbml fields
%   name - name of sbml file
%
% States and parameters consisting of a single character are replaced by
% longer symbols.
% State- and parameter names which coincide with mathematical functions
% in symbolic the Symbolic Toolbox are replaced.
%
% Example:
% arImportSBML('BIOMD0000000379')
%
% Example:
%  ms = arImportSBML('BIOMD0000000379','tend',100)
%  [ms, modelname] = arImportSBML('BIOMD0000000379','tend', 100, 'compartmentbyname')


function varargout = arImportSBMLold(filename, varargin)
warning('This function is deprecated and will be removed in future releases. Use arParseSBML instead.')

%if(~exist('tEnd','var') || isempty(tEnd))
%    tEnd = 100;
%end
%if(~exist('overwrite','var') || isempty(overwrite))
%    overwrite = false;
%end

switches    = {  'tend', 'overwrite', 'keepcompartments', 'compartmentbyname' };
extraArgs   = [       1,           0,                 0,                   0 ];
descriptions = {    { 'Specified tEnd', '' }, ...
                    { 'Overwriting def file', '' }, ...
                    { 'Keeping symbolic compartment names', '' }, ...
                    { 'Using compartment names', '' }, ...
                    };
                
% Parse input arguments
opts = argSwitch( switches, extraArgs, descriptions, 1, varargin );
overwrite = opts.overwrite;

% Set tEnd
tEnd = 100;
if ( opts.tend )
    if ~isnumeric( opts.tend_args )
        error( 'Argument specified after tEnd should be numeric' );
    else
        tEnd = opts.tend_args;
    end
end

if exist('TranslateSBML','file')~=3
    warning('TranslateSBML not found. Please install libSBML and/or add it to Matlab''s search path. EXIT arImportSBML.m now.')
    return
end

% convert relative to absolute path
filename = GetFullPath(filename);

% remove extension
filename = regexprep(filename,'\.(xml|sbml)$','');   

try
    m = TranslateSBML([filename '.xml']);
catch %#ok<CTCH>
    m = TranslateSBML([filename '.sbml']);
end

mIn = m;

m = AdaptVariableNames(m);
m = rules2input(m);
m = findRateRules(m);

% [~,ia] = setdiff({m.species.id},{m.rule.variable});
% m.species = m.species(ia);  % Remove species which are combinations of real dynamic variables, in d2d something like observables

%% check for unsupported features
if(isfield(m,'event') && ~isempty(m.event))
    error('Conversion of events is not yet supported in D2D!')
end
csizes = NaN(size(m.compartment));
for i=1:length(m.compartment)
    csizes(i) = m.compartment(i).size;
    if(sum(m.compartment(i).size<1e-5)>0)
        warning('Small compartment size exist. Could be rescale the equations because of integration problems.')
    end
end

%% model file
[~, name] = fileparts(filename);
new_filename = strrep(name,' ','_');
new_filename = strrep(new_filename,'-','_');
fid = fopen([new_filename '.def'], 'w');

fprintf(fid, 'DESCRIPTION\n');
if(~isempty(m.name))
    fprintf(fid, '"%s"', m.name);
    if(~isempty(m.id))
        fprintf(fid, '" (%s)"', m.id);
    end
    fprintf(fid, '\n');
else
    if(~isempty(m.id))
        fprintf(fid, '"%s"', m.id);
        fprintf(fid, '\n');
    end
end


fprintf(fid, '"SBML level %i version %i"\n', m.SBML_level, m.SBML_version);
if(~isempty(m.notes))
    fprintf(fid, '"%s"\n', m.notes);
end

fprintf(fid, '\nPREDICTOR\n');
if(~isempty(m.time_symbol))
    fprintf(fid, '%s\t T\t "%s"\t time\t 0\t %i\t\n', m.time_symbol, 'n/a', tEnd);
else
    fprintf(fid, 't\t T\t "%s"\t time\t 0\t %i\t\n', 'n/a', tEnd);
end

comps = {};
comp_value = [];
fprintf(fid, '\nCOMPARTMENTS\n');
for j=1:length(m.compartment)
    if(~m.compartment(j).constant)
        error('non-constant compartments are not yet supported in D2D!');
    end
    units = m.compartment(j).units;
    if isempty( units )
        units = 'n/a';
    end
    
    if ( opts.compartmentbyname )
        compName = m.compartment(j).name;
        m.compartmentIDtoD2D = @(j)compartmentIDToName(m,j);
    else
        compName = m.compartment(j).id;
        m.compartmentIDtoD2D = @(j) j;
    end
    if(m.compartment(j).isSetSize)
        if ( opts.keepcompartments )
            fprintf(fid, '%s\t V\t "%s"\t vol.\t \n', sym_check(compName), units);
        else
            fprintf(fid, '%s\t V\t "%s"\t vol.\t %g\n', sym_check(compName), units, m.compartment(j).size);
        end
        comps{end+1} = compName;
        comp_value(end+1) = m.compartment(j).size;
    else
        fprintf(fid, '%s\t V\t "%s"\t vol.\n', sym_check(compName), 'n/a');
    end
end

if ( length(unique(comps)) ~= length(comps) )
    error( 'Duplicate compartment names. Reimport without the flag compartmentbyname.' );
end

fprintf(fid, '\nSTATES\n');
% for j=1:length(m.species)
pat = cell(0); % if length of species names ==1, the species names are extended by '_state'
rep = cell(0);

if isfield(m,'raterule')
    raterulespecies = unique({m.raterule.variable});
    israterule = false(1,length(m.species));
    for i=1:length(raterulespecies)
        israterule = israterule | strcmp(raterulespecies(i),{m.species.id});
    end
else
    israterule = false;
end

for j = find(([m.species.isSetInitialAmount] | [m.species.isSetInitialConcentration]) & ~[m.species.boundaryCondition] | ([m.species.boundaryCondition] & israterule)) % rules should not be defined as states, e.g. K_PP_norm in Huang1996 BIOMD0000000009
    if length(m.species(j).id)==1 %|| strcmp(m.species(j).id,'beta')==1  % special cases or too short
        pat{end+1} =  m.species(j).id; %#ok<AGROW>
        rep{end+1} = [m.species(j).id,'_state']; %#ok<AGROW>
        m.species(j).id2 = rep{end};
        fprintf(fid, '%s\t C\t "%s"\t conc.\t %s\t 1\t "%s"\n', sym_check(rep{end}), 'n/a', ...
            m.compartmentIDtoD2D(m.species(j).compartment), m.species(j).name);
    else  % standard case
        m.species(j).id2 = m.species(j).id;
        fprintf(fid, '%s\t C\t "%s"\t conc.\t %s\t 1\t "%s"\n', sym_check(m.species(j).id), 'n/a', ...
            m.compartmentIDtoD2D(m.species(j).compartment), m.species(j).name);
    end
end

for j=1:length(m.parameter)
    if length(m.parameter(j).id)==1 %|| strcmp(m.species(j).id,'beta')==1  % special cases or too short
        pat{end+1} =  m.parameter(j).id; %#ok<AGROW>
        rep{end+1} = [m.parameter(j).id,'_parameter']; %#ok<AGROW>
        m.parameter(j).id = rep{end};
    end
end

for i=1:length(pat)
    for j=1:length(m.rule)
        m.rule(j).variable = mysubs(m.rule(j).variable,pat{i},rep{i});
        try
            m.rule(j).formula   = mysubs(m.rule(j).formula,pat{i},rep{i});
        catch
            m.rule(j).formula = regexprep(m.rule(j).formula,['(^|(\W)',pat{i},'($|\W)'],['$1',rep{i},'$2'],'all');
        end
    end
end

if(isfield(m,'raterule'))
    for i=1:length(pat)
        for j=1:length(m.raterule)
            m.raterule(j).variable = mysubs(m.raterule(j).variable,pat{i},rep{i});
            try
                m.raterule(j).formula = mysubs(m.raterule(j).formula,pat{i},rep{i});
            catch
                m.raterule(j).formula = regexprep(m.raterule(j).formula,['(^|(\W)',pat{i},'($|\W)'],['$1',rep{i},'$2'],'all');
            end
        end
    end
end

if(isfield(m,'initialAssignment'))
    for i=1:length(pat)
        for j=1:length(m.initialAssignment)
            try
                m.initialAssignment(j).math = mysubs(m.initialAssignment(j).math,pat{i},rep{i});
            catch
                m.initialAssignment(j).math = regexprep(m.initialAssignment(j).math,['(^|(\W)',pat{i},'($|\W)'],['$1',rep{i},'$2'],'all');
            end
        end
    end
end

fprintf(fid, '\nINPUTS\n');
for j=1:length(m.u)
    fprintf(fid, '%s\t C\t "%s"\t conc.\t%s\n', sym_check(m.u(j).variable), m.u(j).units, sym_check(replacePowerFunction(m.u(j).formula)));
end
% treat boundary species as constant inputs
for j=find([m.species.boundaryCondition] & ~israterule)
    m.species(j).id2 = m.species(j).id;
    fprintf(fid, '%s\t C\t "%s"\t conc.\t"%s"\n', sym_check(m.species(j).id), 'n/a', ['init_' m.species(j).id]);
end

arWaitbar(0);
fprintf(fid, '\nREACTIONS\n');
if isfield(m,'raterule')
    for j=1:length(m.raterule)
        arWaitbar(j,length(m.raterule));
        
        prod_spec_name = m.raterule(j).variable;
        for i=1:length(rep)
            prod_spec_name = mysubs(prod_spec_name,pat{i},rep{i});
        end
        prod_spec_name = char(prod_spec_name);
        
        fprintf(fid,'\t -> %s', sym_check(prod_spec_name));
        
        tmpstr = m.raterule(j).formula;
        % repace species names if too short
        for i=1:length(rep)
            tmpstr = mysubs(tmpstr,pat{i},rep{i});
        end
        
        tmpstr = replacePowerFunction(tmpstr);
        
        % replace rules
        
        tmpstr = evalin(symengine,tmpstr);
        findrule = true;
        count = 0;
        while(findrule  && count<100)
            count = count+1;
            
            for jj=1:length(m.rule)
                try
                    tmpstr = mysubs(tmpstr, m.rule(jj).variable, ['(' m.rule(jj).formula ')']);
                catch
                    % rethrow(lasterr)
                    m.rule(jj).variable
                    m.rule(jj).formula
                end
                
            end
            findrule = false;
            vars = symvar(tmpstr);
            for jj=1:length(m.rule)
                if(sum(ismember(vars, sym(m.rule(jj).variable)))>0) %R2013a compatible
                    findrule = true;
                end
            end
        end
        tmpstr = char(tmpstr);
        
        % replace power function
        tmpstr = replacePowerFunction(tmpstr);
        
        fprintf(fid, ' \t CUSTOM "%s" \t"%s"\n', sym_check(tmpstr), m.raterule(j).name);
    end
end

if isfield(m,'reaction') % specified via reactions (standard case)
    for j=1:length(m.reaction)
        arWaitbar(j,length(m.reaction));
        
        % check if reaction constists only of boundary species
        reaction_species = unique({m.reaction(j).reactant(:).species m.reaction(j).product(:).species});
        species_id = NaN(1,length(reaction_species));
        for js=1:length(reaction_species)
            species_id(js) = find(strcmp(reaction_species{js},{m.species.id}));
        end
        if ~all([m.species(species_id).boundaryCondition])
            for jj=1:length(m.reaction(j).reactant)
                % check if reactant is boundary species
                reactant_id = strcmp(m.reaction(j).reactant(jj).species,{m.species.id});
                isboundary = m.species(reactant_id).boundaryCondition;
                if ~isboundary
                    react_spec_name = sym(m.reaction(j).reactant(jj).species);
                    for i=1:length(rep)
                        react_spec_name = mysubs(react_spec_name,pat{i},rep{i});
                    end
                    react_spec_name = char(react_spec_name);
                    
                    if(~isnan(m.reaction(j).reactant(jj).stoichiometry))
                        stoichiometry = m.reaction(j).reactant(jj).stoichiometry;
                    else
                        stoichiometry = 1;
                    end
                    if stoichiometry > 1
                        if(stoichiometry == uint32(stoichiometry))
                            fprintf(fid, '%i %s', stoichiometry, sym_check(react_spec_name));
                        else
                            fprintf(fid, '%2.2f %s', stoichiometry, sym_check(react_spec_name));
                        end
                    else
                        fprintf(fid, '%s', sym_check(react_spec_name));
                    end
                    
                    if(jj ~= length(m.reaction(j).reactant))
                        fprintf(fid, ' + ');
                    end
                    
                end
            end
            if m.reaction(j).reversible
                fprintf(fid, ' \t<-> ');
            else
                fprintf(fid, ' \t-> ');
            end
            for jj=1:length(m.reaction(j).product)
                % check if product is boundary species
                product_id = strcmp(m.reaction(j).product(jj).species,{m.species.id});
                isboundary = m.species(product_id).boundaryCondition;
                if ~isboundary
                    prod_spec_name = sym(m.reaction(j).product(jj).species);
                    for i=1:length(rep)
                        prod_spec_name = mysubs(prod_spec_name,pat{i},rep{i});
                    end
                    prod_spec_name = char(prod_spec_name);
                    
                    if(~isnan(m.reaction(j).product(jj).stoichiometry))
                        stoichiometry = m.reaction(j).product(jj).stoichiometry;
                    else
                        stoichiometry = 1;
                    end
                    if(stoichiometry~=1)
                        if(stoichiometry==uint32(stoichiometry))
                            fprintf(fid, '%i %s', stoichiometry, sym_check(prod_spec_name));
                        else
                            fprintf(fid, '%2.2f %s', stoichiometry, sym_check(prod_spec_name));
                        end
                    else
                        fprintf(fid, '%s', sym_check(prod_spec_name));
                    end
                    if(jj ~= length(m.reaction(j).product))
                        fprintf(fid, ' + ');
                    end
                    
                end
            end
            
            tmpstr = arSym(m.reaction(j).kineticLaw.math);
            % repace species names if too short
            for i=1:length(rep)
                tmpstr = armysubs(tmpstr,pat{i},rep{i});
            end
            
            % make parameters unique
            if(isfield(m.reaction(j).kineticLaw, 'parameter'))
                for jj=1:length(m.reaction(j).kineticLaw.parameter)
                    tmpstr = mysubs(tmpstr, m.reaction(j).kineticLaw.parameter(jj).id, ...
                        [m.reaction(j).id '_' m.reaction(j).kineticLaw.parameter(jj).id]);
                end
            end
            
            % divide rates by compartment volume
            reaction_comp = findReactionCompartment(m,j, csizes);
            reaction_comp = m.compartmentIDtoD2D( reaction_comp );
            
            if ~isempty(reaction_comp) %&& sum(strcmp(reaction_comp,strsplit(char(tmpstr),'*')))==1
                tmpstr = ['(' char(tmpstr) ')/' sym_check(reaction_comp)];
            end
            
            % replace compartment volumes if requested
            if (~opts.keepcompartments)
                for jj=1:length(m.compartment)
                    tmpstr = mysubs(tmpstr, m.compartmentIDtoD2D( m.compartment(jj).id ), num2str(m.compartment(jj).size));
                end
            else
                for jj=1:length(m.compartment)
                    tmpstr = mysubs(tmpstr, m.compartmentIDtoD2D( m.compartment(jj).id ), sprintf('vol_%s', m.compartmentIDtoD2D( m.compartment(jj).id ) ));
                end
            end
            
            % remove functions
            tmpstr = char(tmpstr);
            for jj=1:length(m.functionDefinition)
                tmpfun = m.functionDefinition(jj).math;
                tmpfun = strrep(tmpfun, 'lambda(', '');
                tmpfun = tmpfun(1:end-1);
                tmpfun = replacePowerFunction(tmpfun,false);
                C = textscan(tmpfun, '%s', 'Whitespace', ',');
                C = C{1};
                
                tmpstr = replaceFunction(tmpstr, m.functionDefinition(jj).id, C(1:end-1), C(end));
                
            end
            
            % replace power function
            tmpstr = replacePowerFunction(tmpstr);
            
            % replace rules
            tmpstr = sym(tmpstr);
            findrule = true;
            count = 0;
            while(findrule && count < 100)
                count = count+1;
                for jj=1:length(m.rule)
                    tmpstr = mysubs(tmpstr, m.rule(jj).variable, ['(' m.rule(jj).formula ')']);
                end
                findrule = false;
                vars = symvar(tmpstr);
                for jj=1:length(m.rule)
                    if(sum(ismember(vars, sym(m.rule(jj).variable)))>0) %R2013a compatible
                        findrule = true;
                    end
                end
            end
            tmpstr = char(tmpstr);
            
            tmpstr = replacePowerFunction(tmpstr);
            
            fprintf(fid, ' \t CUSTOM "%s" \t"%s"\n', sym_check(tmpstr), m.reaction(j).id);
        end
    end
end   % end if dynamics specified either raterule or reaction

arWaitbar(-1);

fprintf(fid, '\nDERIVED\n');

fprintf(fid, '\nCONDITIONS\n');
% compartments
if ( opts.keepcompartments )
    for j = 1 : length( m.compartment )
        fprintf( fid, 'vol_%s   "%g"\n', m.compartmentIDtoD2D( m.compartment(j).id ), m.compartment(j).size );
    end
end

fprintf(fid, '\nPARAMETERS\n');
specs = {};
spec_value = [];
for j=1:length(m.species)
    if(m.species(j).isSetInitialConcentration)
        ub = 1000;
        if(m.species(j).initialConcentration>ub)
            ub = m.species(j).initialConcentration*10;
        end
        fprintf(fid, 'init_%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.species(j).id2), ...
            m.species(j).initialConcentration, m.species(j).constant==0, ub);
        specs{end+1} = m.species(j).id2;
        spec_value(end+1) = m.species(j).initialConcentration;
    elseif(m.species(j).isSetInitialAmount)
        comp_id = strcmp(m.species(1).compartment,{m.compartment.id});
        comp_vol = m.compartment(comp_id).size;
        initial_conc = m.species(j).initialAmount/comp_vol;
        ub = 1000;
        if(initial_conc>ub)
            ub = initial_conc*10;
        end
        fprintf(fid, 'init_%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.species(j).id2), ...
            initial_conc, m.species(j).constant==0, ub);
        specs{end+1} = m.species(j).id2;
        spec_value(end+1) = initial_conc;
    end
end
pars = {};
par_value = [];
for j=1:length(m.parameter)
    isRule = false;
    for jjj=1:length(m.rule)
        if(strcmp(m.rule(jjj).variable, m.parameter(j).id))
            isRule = true;
        end
        if(strcmp(m.rule(jjj).variable, m.parameter(j).name))
            isRule = true;
        end
    end
    if(~isRule)
        if(m.parameter(j).isSetValue)
            ub = 1000;
            if(m.parameter(j).value>ub)
                ub = m.parameter(j).value*10;
            end
            fprintf(fid, '%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.parameter(j).id), ...
                m.parameter(j).value, m.parameter(j).constant==0, ub);
            pars{end+1} = m.parameter(j).id;
            par_value(end+1) = m.parameter(j).value;
            
        end
    end
end
% parameters from reactions
for j=1:length(m.reaction)
    if(isfield(m.reaction(j).kineticLaw, 'parameter'))
        for jj=1:length(m.reaction(j).kineticLaw.parameter)
            if(m.reaction(j).kineticLaw.parameter(jj).isSetValue)
                ub = 1000;
                if(m.reaction(j).kineticLaw.parameter(jj).value>ub)
                    ub = m.reaction(j).kineticLaw.parameter(jj).value*10;
                end
                fprintf(fid, '%s_%s\t %g\t %i\t 0\t 0\t %g\n', m.reaction(j).id, ...
                    m.reaction(j).kineticLaw.parameter(jj).id, ...
                    m.reaction(j).kineticLaw.parameter(jj).value, m.reaction(j).kineticLaw.parameter(jj).constant==0, ub);
            end
        end
    end
end

% initial assignments
for i=1:length(m.initialAssignment)
    assignment_value = m.initialAssignment(i).math;
    assignment_value = subs(assignment_value, pars, par_value);
    assignment_value = subs(assignment_value, specs, spec_value);
    assignment_value = subs(assignment_value, comps, comp_value);
    assignment_value = eval(assignment_value);
    
    if any(strcmp(m.initialAssignment(i).symbol,{m.parameter.id}))
        ub = 1000;
        if assignment_value > ub
            ub = 10*assignment_value;
        end
        fprintf(fid, '%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.initialAssignment(i).symbol), ...
            assignment_value, 1, ub);
    end
end

fclose(fid);

%% data file

fid = fopen([new_filename '_data.def'], 'w');

fprintf(fid, 'DESCRIPTION\n');
if(~isempty(m.name))
    fprintf(fid, '"data file for %s"\n', m.name);
end

fprintf(fid, '\nPREDICTOR\n');
if(~isempty(m.time_symbol))
    fprintf(fid, '%s\t T\t "%s"\t time\t 0\t %i\t\n', m.time_symbol, 'n/a', 100);
else
    fprintf(fid, 't\t T\t "%s"\t time\t 0\t %i\t\n', 'n/a', 100);
end

fprintf(fid, '\nINPUTS\n');

% TODO: Access parameter name list. Where to implement and how?
% TODO: Access Units from model. 

if isfield(m,'rule') && ~isempty(m.rule(1).formula) % look for observation functions
    for j = 1:length(m.rule)
        split_err = strsplit(m.rule(j).variable,'sigma_');
        split_obs = strsplit(m.rule(j).variable,'observable_');
        isobs(j) = length(split_obs) == 2;
        iserr(j) = length(split_err) == 2;
        if isobs(j)
            obsname{j} = split_obs{2};
        elseif iserr(j)
            obsname{j} = split_err{2};
        else
            obsname{j} = '';
        end
    end
    [~,B,C] = unique(obsname);
    jsobs = find(isobs);
    jserr = find(iserr);
    
    fprintf(fid, '\nOBSERVABLES\n');
    for j=1:length(B)
        idx = find((j == C) & isobs');
        if ~isempty(idx)
        fprintf(fid, '%s \t C\t "%s"\t "conc."\t 0 0 "%s" "%s"\n', sym_check(m.rule(idx).variable), 'n/a', ...
            m.rule(idx).formula, m.rule(idx).name); 
        end
    end
    
    fprintf(fid, '\nERRORS\n');
    for j=1:length(B)
        idxobs = find((j == C) & isobs');
        idxerr = find((j == C) & iserr');
        if (length(idxobs) == 1) && (length(idxerr) == 1)
            fprintf(fid, '%s\t "%s"\n', sym_check(m.rule(idxobs).variable), sym_check(m.rule(idxerr).formula));
        else
            warning('Double check error model functions in data.def. Something might have gone wrong.')
        end
    end
elseif isfield(m.species,'id2')%old behavior
    fprintf(fid, '\nOBSERVABLES\n');
    for j=1:length(m.species)
        fprintf(fid, '%s_obs\t C\t "%s"\t conc.\t 0 0 "%s" "%s"\n', sym_check(m.species(j).id2), 'n/a', ...
            m.species(j).id2, m.species(j).name);
    end
    
    fprintf(fid, '\nERRORS\n');
    for j=1:length(m.species)
        fprintf(fid, '%s_obs\t "sd_%s"\n', sym_check(m.species(j).id2), sym_check(m.species(j).id2));
    end
    warning('No real observation functions defined!')
end

fprintf(fid, '\nCONDITIONS\n');

fprintf(fid, '\nPARAMETERS\n');

fclose(fid);

if ~isdir('./Models')
    mkdir('Models');
end
if(overwrite)
    movefile([new_filename '.def'],'Models','f');
%     system(['mv -f ',new_filename '.def Models']);
else
    dest = ['Models',filesep,new_filename '.def'];
    if exist(dest,'file')==0
        movefile([new_filename '.def'],'Models');
    else
        fprintf('%s already exists. Either use the flag ''overwrite'' or move the files by hand.\n',dest);
    end
%     system(['mv ',new_filename '.def Models']);
end

if ~isdir('./Data')
    mkdir('Data');
end
if(overwrite)    
    movefile([new_filename '_data.def'],'Data','f');
%     system(['mv -f ',new_filename '_data.def Data']);
else
    dest = ['Data',filesep,new_filename '_data.def'];
    if exist(dest,'file')==0
        movefile([new_filename '_data.def'],'Data');
    else
        fprintf('%s already exists. Either use the flag ''overwrite'' or move the files by hand.\n',dest);
    end
%     system(['mv ',new_filename '_data.def Data']);
end

% generate Setup.m
if(~exist('Setup.m','file'))
    fprintf('Generate a standard template for loading and compiling the model (Setup.m) ...\n')
    fid = fopen('Setup.m','w');
    fprintf(fid, 'arInit;\n');
    fprintf(fid, 'arLoadModel(''%s'');\n',new_filename);
    fprintf(fid, 'arLoadData(''%s'');\n',[new_filename '_data']);
    fprintf(fid, 'arCompileAll;\n');
else
    fprintf('Setup.m already available in the working directory.\n')
end

fprintf('Model- and data definition files created.\n');
fprintf('If loading via Setup.m does not work, try to solve issues by adapting the def files by hand.\n');
fprintf('The following issues might occur:\n')
fprintf(' - In SBML, compartements are often used as a kind of annotation and not for considering volume factors. Then compartement (should) have the same volumes. In D2D, a reaction of compounds located in different compartements might cause an error.\n');
fprintf(' - Non-standard function, e.g. something like "stepfunction()", which are not available in C have to be replaced or implemented properly.\n')
fprintf(' - The time axis (termed "PREDICTOR") might not assigned properly. Replace the respective variables properly.\n')

% assign returns
if nargout > 0
    ms.d2d = m;
    ms.sbml = mIn;
    ms.pat = pat;
    ms.rep = rep;
    varargout{1} = ms;
end
if nargout > 1
    varargout{2} = new_filename;
end


function str = replaceFunction(str, funstr, C, funmat)
% %% test replace functions
% clc
% str = 'k1 + power(k1*2, k2+(7*log(k3))) + 10*p3 + power(k1*2, k2+(7*log(k3))) + 10*p3';
% funstr = 'power';
% C{1} = 'a';
% C{2} = 'b';
% funmat = 'a^b';
% disp(str);

funindex = strfind(str, [funstr '(']);
while(~isempty(funindex))
    
    substr = str(funindex(1):end);
    
    openindex = strfind(substr, '(');
    closeindex = strfind(substr, ')');
    
    mergedindex = [openindex closeindex];
    rankingindex = [ones(size(openindex)) -ones(size(closeindex))];
    
    [sortedmergedindex, isortedindex] = sort(mergedindex);
    sortedrankingindex = rankingindex(isortedindex);
    
    endfunindex = find(cumsum(sortedrankingindex)==0);
    if(isempty(endfunindex))
        error('bracketing error close to function %s', funstr);
    end
    endfunindex = sortedmergedindex(endfunindex(1));
    
    substr = substr(openindex+1:endfunindex-1);
    
    D = textscan(substr, '%s', 'Whitespace', ',');
    D = D{1};
    if(length(C)~=length(D))
        error('input output parameter mismatch');
    end
    
    funtmplate = funmat;
    
    % Replace longest names first               %#<JV>
    [~,I]=sort(cellfun(@length,C), 'descend');  %#<JV>
    for j=1:length(D)
        pattern = sprintf('(^%s|%s$)|((?<=[\\(\\+\\*\\-\\/])(%s)(?=[\\)\\+\\*\\-\\/]))', C{I(j)}, C{I(j)}, C{I(j)}); % use regex for replacing pars and vars correctly
        funtmplate = regexprep(funtmplate, pattern, ['(' D{I(j)} ')']);
    end
    funtmplate = ['(' funtmplate ')']; %#ok<AGROW>
    % disp(funtmplate)
    
    if(funindex(1)-1>1 && funindex(1)+endfunindex<length(str))
        str = [str(1:funindex(1)-1) funtmplate str(funindex(1)+endfunindex:end)];
    elseif(funindex(1)-1>1)
        str = [str(1:funindex(1)-1) funtmplate];
    elseif(funindex(1)+endfunindex<length(str))
        str = [funtmplate str(funindex(1)+endfunindex:end)];
    else
        str = funtmplate;
    end
    str = cell2mat(str);
    % disp(str)
    
    funindex = strfind(str, funstr);
end
% disp(str)

str = char(sym(str));


function str = replacePowerFunction(str, issym)
% replace power function ('power' and 'pow')
% str = 'k1 + power(k1*2, k2+(7*log(k3))) + 10*p3 + power(k1*2, k2+(7*log(k3))) + 10*p3';
%
% issym:    set to true, when used togehter with symbolic evaluation.
% Replaces 'power' with '_power'

narginchk(1,2)
if(~exist('issym','var'))
    issym = true;
end

if issym
    str = strrep(str, 'power(', '_power('); % FIXME: use regexp instead
    str = strrep(str, 'pow(', '_power('); % FIXME: use regexp instead
else
    C = {'a','b'};
    funstr = 'power';
    
    str = char(str);
    % disp(str);
    funindex = strfind(str, [funstr '(']);
    while(~isempty(funindex))
        
        substr = str(funindex(1):end);
        
        openindex = strfind(substr, '(');
        closeindex = strfind(substr, ')');
        
        mergedindex = [openindex closeindex];
        rankingindex = [ones(size(openindex)) -ones(size(closeindex))];
        
        [sortedmergedindex, isortedindex] = sort(mergedindex);
        sortedrankingindex = rankingindex(isortedindex);
        
        endfunindex = find(cumsum(sortedrankingindex)==0);
        if(isempty(endfunindex))
            error('bracketing error close to function %s', funstr);
        end
        endfunindex = sortedmergedindex(endfunindex(1));
        
        substr = substr(openindex+1:endfunindex-1);
        
        D = textscan(substr, '%s', 'Whitespace', ',');
        D = D{1};
        if(length(C)~=length(D))
            error('input output parameter mismatch');
        end
        
        funtmplate = sprintf('((%s)^(%s))',D{1},D{2});
        %     disp(funtmplate)
        
        if(funindex(1)-1>1 && funindex(1)-1+endfunindex<length(str)) % in between
            str = [str(1:funindex(1)-1) funtmplate str(funindex(1)+endfunindex:end)];
        elseif(funindex(1)-1>1) % at begining
            str = [str(1:funindex(1)-1) funtmplate];
        elseif(funindex(1)-1+endfunindex<length(str)) % at end
            str = [funtmplate str(funindex(1)+endfunindex:end)];
        else % whole string
            str = funtmplate;
        end
        %     disp(str)
        
        funindex = strfind(str, funstr);
    end
end


function m = findRateRules(m)

drin = [];
for i=1:length(m.rule)
    switch m.rule(i).typecode
        case 'SBML_ASSIGNMENT_RULE'
            drin = [drin,i];%#ok<AGROW> % standard case
        case 'SBML_RATE_RULE'
            if ~isfield(m,'raterule')
                m.raterule = m.rule(i);
            else
                m.raterule(end+1) = m.rule(i);
            end
        otherwise
            m.rule(i).typecode
            error(' m.rule(i).typecode unknown');
    end
    
end
m.rule = m.rule(drin);


function m = rules2input(m)
is_input = zeros(size(m.rule));
for r=1:length(m.rule)
    s = symvar(m.rule(r).formula);
    if(~isempty(intersect(s,'TIME')))
        is_input(r) = 1;
    end
end

m.u = m.rule(is_input==1);
m.rule = m.rule(is_input~=1);

for i=1:length(m.u)
    %     m.u(i).formula = char(mysubs(sym(m.u(i).formula),'time','t')); % does
    %     not work, at least in R2014a
    m.u(i).formula = strrep(m.u(i).formula,'TIME','t');
end


function m = AdaptVariableNames(m)
% The following function will
%   alter variable names which cannot be handled by the Symbolic Math
%   function subs()

for i=1:length(m.species)
    m.species(i).id = sym_check(m.species(i).id);
    m.species(i).compartment = sym_check(m.species(i).compartment);
end


for i=1:length(m.parameter)
    m.parameter(i).id = sym_check(m.parameter(i).id);
end

for i=1:length(m.rule)
    m.rule(i).variable = sym_check(m.rule(i).variable);
    m.rule(i).formula = sym_check(m.rule(i).formula);
end

for i=1:length(m.reaction)
    for j=1:length(m.reaction(i).reactant)
        m.reaction(i).reactant(j).species = sym_check(m.reaction(i).reactant(j).species);
    end
    m.reaction(i).kineticLaw.math = sym_check( m.reaction(i).kineticLaw.math);
    if isfield(m.reaction(i).kineticLaw,'formula')
        m.reaction(i).kineticLaw.formula = sym_check( m.reaction(i).kineticLaw.formula);
    end
end


function s = sym_check(s)
% Replacement if symbolic variable coincides with function. Without
% replacement subs would not work.
keywords = {'time','gamma','sin','cos','tan','beta','log','asin','atan','acos','acot','cot','theta','D','I','E'};

issym = strcmp(class(s),'sym'); %#ok<STISA>
if(issym)
    s = char(s);
end

sv = symvar(s);

svinter = intersect(sv,keywords);

for i=1:length(svinter)
    for j=1:3
        if(length(svinter{i})>1)
            s = regexprep(s,['(^|[\W])',svinter{i},'([\W]|$)'],['$1',upper(svinter{i}),'$2'],'all');
        else
            s = regexprep(s,['(^|[\W])',svinter{i},'([\W]|$)'],['$1',svinter{i},'_symbol','$2'],'all');
        end
    end
end

if(issym)
    s = sym(s);
end


function s = mysubs(s,pat,rep)

keywords = {'time','gamma','sin','cos','tan','beta','log','asin','atan','acos','acot','cot','theta','D','I','E'};

issym = strcmp(class(s),'sym'); %#ok<STISA>
if(~issym)
    s = evalin(symengine,s);
end

issym = strcmp(class(pat),'sym'); %#ok<STISA>
if(~issym)
    pat_sym = evalin(symengine,pat);
else
    pat_sym = pat;
end

issym = strcmp(class(rep),'sym'); %#ok<STISA>
if(~issym)
    rep_sym = evalin(symengine,rep);
else
    rep_sym = rep;
end

if isempty(intersect(pat,keywords))
    try
        s = char(subs(s,pat_sym,rep_sym));
        err=0;
    catch lasterr
        disp(lasterr)
        err=1;
    end
end

if ~isempty(intersect(pat,keywords)) || err==1
    % symbolic toolbox keywords (function) do not work with subs
    sv = symvar(char(s));
    svcell = cell(size(sv));
    for i=1:length(sv)
        svcell{i} = char(sv(i));
    end
    if ~isempty(intersect(pat,svcell))
        s = char(s);
        for j=1:3
            s = regexprep(s,['(^|(\W)',pat,'($|\W)'],['$1',rep,'$2'],'all');
        end
        if(issym)
            s = sym(s);
        end
    end
end

function c = compartmentIDToName( m, reaction_comp )
    for jc = 1 : length( m.compartment )
        if ( strcmp( reaction_comp, m.compartment(jc).id ) )
            c = m.compartment(jc).name;
            return;
        end
    end

function c = findReactionCompartment(m, j, csizes)
% find compartment of reacting species to convert from SBML rate convention
% (particle flux) to d2d (concentration flux).

comp_r = {};
for jr=1:length(m.reaction(j).reactant);
    js = strcmp({m.species.id},m.reaction(j).reactant(jr).species);
    comp_r{jr} = m.species(js).compartment; %#ok<AGROW>
end

comp_p = {};
for jr=1:length(m.reaction(j).product);
    js = strcmp({m.species.id},m.reaction(j).product(jr).species);
    comp_p{jr} = m.species(js).compartment; %#ok<AGROW>
end

% check educt and product compartements for consistency
if ~isempty(comp_r)
    if length(unique(comp_r))~=1
        if length(unique(csizes))==1
            warning('Reactants originate from more than one compartment. Such a model definition is in general not reasonable but does not matter in this case because all compartments have the same size.');
        else
            error('Reactants originate from more than one compartment. Such a model definition is in general not reasonable.');            
        end
    end
end
if ~isempty(comp_p)
    if length(unique(comp_p))~=1
        if length(unique(csizes))==1
            warning('Products are in more than one compartment. Such a model definition is in general not reasonable but does not matter in this case because all compartments have the same size.');
        else
            error('Products are in more than one compartment. Such a model definition is not reasonable.');
        end
    end
end

% educt and product exist. Have to divide by source volume to get D2D
% convention.
% c = comp_r{1};

if isempty(comp_r) && ~isempty(comp_p)
     c = comp_p{1};
else
    c = comp_r{1};
end
% 
% if ~isempty(comp_r) && ~isempty(comp_p)
%     % educt and product in the same compartment
%     if isequal(unique(comp_r),unique(comp_p))
%         c = comp_r{1};
%     end
%     % only educt has a compartment
% elseif ~isempty(comp_r) && isempty(comp_p)
%     c = comp_r{1};
%     % only product has a compartment
% elseif isempty(comp_r) && ~isempty(comp_p)
%     c = comp_p{1};
% end

