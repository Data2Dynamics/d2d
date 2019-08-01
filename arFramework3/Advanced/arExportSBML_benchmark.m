% function arExportSBML_benchmark(m, d, steadystate)
%
% export model to SBML using libSBML
%
% m:            model index
% d:            data index
% steadystate:  prequilibrate condition as in ar file
%
%%% Differences to arExportSBML: %%%
% Observation functions are written out via assignment rules. If the
% observation function is log-transformed in ar, the log trafo will be
% written into SBML as well.

function arExportSBML_benchmark(m, d, steadystate)

global ar

c = ar.model(m).data(d).cLink;
% simulate once for initial values
arSimu(0,1,0)

copasi = true;

if(~exist('data','var'))
    data = 1;
end

if (~exist('steadystate','var'))
    if ( isfield( ar.model(m), 'ss_condition' ) )
        steadystate = true;
    else
        steadystate = false;
    end
end

try
M = TranslateSBML(which('empty.xml'));
F = TranslateSBML(which('filled.xml'));
catch
    warning('error in libSBML. Probably backwards compatibility issues with old MATLAB version. Should work with 2019a')
end

M.id = ar.model(m).name;
if(~isempty(ar.model(m).description))
    M.notes = ar.model(m).description{1};
end
%% compartements
if(~isempty(ar.model(m).c))
    for jc = 1:length(ar.model(m).c)
        M.compartment(jc).typecode = 'SBML_COMPARTMENT';
        M.compartment(jc).metaid = '';
        M.compartment(jc).notes = '';
        M.compartment(jc).annotation = '';
        M.compartment(jc).sboTerm = -1;
        M.compartment(jc).name = '';
        M.compartment(jc).id = ar.model(m).c{jc};
        M.compartment(jc).compartmentType = '';
        M.compartment(jc).spatialDimensions = 3;
        M.compartment(jc).constant = 1;
        if ~isempty(ar.model(m).cUnits)
            M.compartment(jc).units = ar.model(m).cUnits{jc,2};
        else
            M.compartment(jc).units = '';
        end
        M.compartment(jc).outside = '';
        M.compartment(jc).isSetSize = 1;
        M.compartment(jc).isSetVolume = 1;
        M.compartment(jc).level = 2;
        M.compartment(jc).version = 4;
        
        qp = ismember(ar.pLabel, ar.model(m).pc{jc}); %R2013a compatible
        if(sum(qp)==1)
            pvalue = ar.p(qp);
            if(ar.qLog10(qp))
                pvalue = 10^pvalue;
            end
            M.compartment(jc).size = pvalue;
        elseif(sum(qp)==0)
            qp = ismember(ar.model(m).condition(c).pold, ar.model(m).pc{jc}); %R2013a compatible
            if(sum(qp)==1)
                pvalue = ar.model(m).condition(c).fp{qp};
                M.compartment(jc).size = str2num(pvalue); %#ok str2num also identifies brackets, in which case str2double simply returns NaN
            else
                pvalue = str2num(ar.model(m).pc{jc}); %#ok
                if(~isnan(pvalue))
                    M.compartment(jc).size = pvalue;
                else
                    error('%s not found', ar.model(m).pc{jc});
                end
            end
        else
            error('%s not found', ar.model(m).pc{jc});
        end
    end
else
    M.compartment(1).typecode = 'SBML_COMPARTMENT';
    M.compartment(1).metaid = '';
    M.compartment(1).notes = '';
    M.compartment(1).annotation = '';
    M.compartment(1).sboTerm = -1;
    M.compartment(1).name = '';
    M.compartment(1).id = 'default';
    M.compartment(1).compartmentType = '';
    M.compartment(1).spatialDimensions = 3;
    M.compartment(1).constant = 1;
    if ~isempty(ar.model(m).cUnits)
        M.compartment(1).units = ar.model(m).cUnits{1,2};
    else
        M.compartment(1).units = '';
    end
    M.compartment(1).outside = '';
    M.compartment(1).isSetSize = 1;
    M.compartment(1).isSetVolume = 1;
    M.compartment(1).level = 2;
    M.compartment(1).version = 4;
    M.compartment(1).size = 1;
end

%% species
Crules = {};

simulated_ss = 0;
if ( steadystate )
    if (isfield(ar.model(m).condition(c),'ssLink') && ~isempty( ar.model(m).condition(c).ssLink ) )
        warning( 'Using simulated steady state values as initial condition (non-parametric)' );
        x_ss = ar.model(m).ss_condition(ar.model(m).condition(c).ssLink).xFineSimu(end,:);
        simulated_ss = 1;
    end
end
for jx = 1:length(ar.model(m).x)
    M.species(jx).typecode = 'SBML_SPECIES';
    M.species(jx).metaid = '';
    M.species(jx).notes = '';
    M.species(jx).annotation = '';
    M.species(jx).sboTerm = -1;
    if ( isfield( ar.model(m), 'xNames' ) )
        M.species(jx).name = ar.model(m).xNames{jx};
    else
        M.species(jx).name = '';
    end
    M.species(jx).id = ar.model(m).x{jx};
    M.species(jx).speciesType = '';
    if(length(ar.model(m).cLink)>=jx && ~isempty(ar.model(m).cLink(jx)) && ar.model(m).cLink(jx)~=0)
        M.species(jx).compartment = ar.model(m).c{ar.model(m).cLink(jx)};
    else
        M.species(jx).compartment = 'default';
    end
    M.species(jx).initialAmount = NaN;
    if ~isempty(ar.model(m).xUnits)
        M.species(jx).substanceUnits = ar.model(m).xUnits{jx,2};
    else
        M.species(jx).substanceUnits = '';
    end
    M.species(jx).hasOnlySubstanceUnits = 0;
    M.species(jx).boundaryCondition = 0;
    M.species(jx).charge = 0;
    M.species(jx).constant = 0;
    M.species(jx).isSetInitialAmount = 0;
    M.species(jx).isSetInitialConcentration = 1;
    M.species(jx).isSetCharge = 0;
    M.species(jx).level = 2;
    M.species(jx).version = 4;
    
    qp = ismember(ar.pLabel, ar.model(m).px0{jx}); %R2013a compatible

    if ( ~simulated_ss )
        % check if init parameter still exists in condition parameters
        is_set_cond = sum(ismember(ar.model(m).condition(c).fp, ar.model(m).px0{jx}))==0;
        if(sum(qp)==1 && ~is_set_cond)
            M.species(jx).initialConcentration = 1;
            Crules{end+1,1} = ar.model(m).x{jx}; %#ok<AGROW>
            Crules{end,2} = ar.pLabel{qp}; %#ok<AGROW>
        elseif(sum(qp)==0 || is_set_cond)
            qp = ismember(ar.model(m).condition(c).pold, ar.model(m).px0{jx}); %R2013a compatible
            if(sum(qp)==1)
                pvalue = char(arSym(ar.model(m).condition(c).fp{qp}));
                if(~isnan(str2num(pvalue))) %#ok
                    pvalue = str2num(pvalue); %#ok
                    M.species(jx).initialConcentration = pvalue;
                else
                    Crules{end+1,1} = ar.model(m).x{jx}; %#ok<AGROW>
                    Crules{end,2} = pvalue; %#ok<AGROW>
                    M.species(jx).initialConcentration = 1;
                end
            else
                error('%s not found', ar.model(m).pc{jc});
            end
        else
            error('%s not found', ar.model(m).pc{jc});
        end
    else
        M.species(jx).initialConcentration = x_ss(jx);
    end
end

%% parameters
% for jp = 1:length(ar.model(m).condition(c).p)
%     M.parameter(jp).typecode = 'SBML_PARAMETER';
%     M.parameter(jp).metaid = '';
%     M.parameter(jp).notes = '';
%     M.parameter(jp).annotation = '';
%     M.parameter(jp).sboTerm = -1;
%     M.parameter(jp).name = '';
%     M.parameter(jp).id = ar.model(m).condition(c).p{jp};
%     M.parameter(jp).units = '';
%     M.parameter(jp).constant = 1;
%     M.parameter(jp).isSetValue = 1;
%     M.parameter(jp).level = 2;
%     M.parameter(jp).version = 4;
%     
%     qp = ismember(ar.pLabel, ar.model(m).condition(c).p{jp}); %R2013a compatible
%     if(sum(qp)==1)
%         pvalue = ar.p(qp);
%         if(ar.qLog10(qp) == 1)
%             pvalue = 10^pvalue;
%         end
%     else
%         pvalue = 1;
%     end
%     M.parameter(jp).value = pvalue;
% end

for id = 1:length(ar.model(m).data(d).p)
    %id_tmp = id+length(ar.model(m).condition(c).p);
    id_tmp = id;
    M.parameter(id_tmp).typecode = 'SBML_PARAMETER';
    M.parameter(id_tmp).metaid = '';
    M.parameter(id_tmp).notes = '';
    M.parameter(id_tmp).annotation = '';
    M.parameter(id_tmp).sboTerm = -1;
    M.parameter(id_tmp).name = ar.model(m).data(d).p{id};                       
    M.parameter(id_tmp).id = ar.model(m).data(d).p{id};
    M.parameter(id_tmp).units = '';
    M.parameter(id_tmp).constant = 0;
    M.parameter(id_tmp).isSetValue = 1;
    M.parameter(id_tmp).level = 2;
    M.parameter(id_tmp).version = 4;
    
    qp = ismember(ar.pLabel, ar.model(m).data(d).p{id}); %R2013a compatible
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp) == 1)
            pvalue = 10^pvalue;
        end
    else
        pvalue = 1;
    end
    M.parameter(id_tmp).value = pvalue;
end

for iY = 1:length(ar.model(m).data(d).y)
    id_tmp = length(M.parameter)+1;
    M.parameter(id_tmp).typecode = 'SBML_PARAMETER';
    M.parameter(id_tmp).metaid = '';
    M.parameter(id_tmp).notes = '';
    M.parameter(id_tmp).annotation = '';
    M.parameter(id_tmp).sboTerm = -1;
    M.parameter(id_tmp).name = '';
    M.parameter(id_tmp).id = ar.model(m).data(d).y{iY};
    M.parameter(id_tmp).units = '';
    M.parameter(id_tmp).constant = 0;
    M.parameter(id_tmp).isSetValue = 1;
    M.parameter(id_tmp).level = 2;
    M.parameter(id_tmp).version = 4;    
    M.parameter(id_tmp).value = 0;
end

%% rules
for jr = 1:size(Crules,1)
    M.initialAssignment(jr).typecode = 'SBML_INITIAL_ASSIGNMENT';
    M.initialAssignment(jr).metaid = '';
    M.initialAssignment(jr).notes = '';
    M.initialAssignment(jr).annotation = '';
    M.initialAssignment(jr).sboTerm = -1;
    M.initialAssignment(jr).symbol = Crules{jr,1};
    M.initialAssignment(jr).math = Crules{jr,2};
    M.initialAssignment(jr).level = 2;
    M.initialAssignment(jr).version = 4;
end

%% reactions
fv = ar.model(m).fv;
fv = arSym(fv);
% fv = arSubs(fv, ar.model(m).u, ar.model(m).condition(c).fu');
fv = arSubs(fv, arSym(ar.model(m).condition(c).pold), arSym(ar.model(m).condition(c).fp'));
fv = arSubs(fv, arSym(ar.model(m).t), arSym('time'));

vcount = 1;
arWaitbar(0);
for jv = 1:length(ar.model(m).fv)
    arWaitbar(jv, length(ar.model(m).fv));
    ratetemplate = fv(jv);
    
    if(ratetemplate~=0)
        M.reaction(vcount).typecode = 'SBML_REACTION';
        M.reaction(vcount).metaid = '';
        M.reaction(vcount).notes = '';
        M.reaction(vcount).annotation = '';
        M.reaction(vcount).sboTerm = -1;
        if(isfield(ar.model(m),'v') && length(ar.model(m).v)>=jv && ~isempty(ar.model(m).v))
            M.reaction(vcount).name = ar.model(m).v{jv};
        else
            M.reaction(vcount).name = '';
        end
        if ( isfield( ar.model(m), 'reversible' ) && ~isempty(ar.model(m).reversible))
            M.reaction(vcount).reversible = ar.model(m).reversible(jv);
        else
            M.reaction(vcount).reversible = 0;
        end
        M.reaction(vcount).fast = 0;%-1;
        M.reaction(vcount).isSetFast = 0;
        M.reaction(vcount).level = 2;
        M.reaction(vcount).version = 4;
        
        if(isfield(ar.model(m),'v') && ~isempty(ar.model(m).v))
            % replace spaces with underscores
            M.reaction(vcount).id = sprintf( 'v%d_%s', jv, strrep(ar.model(m).v{jv},' ','_') );
        else
            M.reaction(vcount).id = sprintf('reaction%i', jv);
        end
        
        %set empty struct for reactant
        M.reaction(vcount).reactant = F.reaction(1).reactant;
        scount = 1;
        scomp = [];
        for jsource = find(ar.model(m).N(:,jv)<0)'
            M.reaction(vcount).reactant(scount).typecode = 'SBML_SPECIES_REFERENCE';
            M.reaction(vcount).reactant(scount).metaid = '';
            M.reaction(vcount).reactant(scount).notes = '';
            M.reaction(vcount).reactant(scount).annotation = '';
            M.reaction(vcount).reactant(scount).sboTerm = -1;
            M.reaction(vcount).reactant(scount).species = ar.model(m).x{jsource};
            M.reaction(vcount).reactant(scount).id = '';
            M.reaction(vcount).reactant(scount).name = '';
            M.reaction(vcount).reactant(scount).stoichiometry = abs(ar.model(m).N(jsource,jv));
            M.reaction(vcount).reactant(scount).stoichiometryMath = F.reaction(2).reactant.stoichiometryMath;
            M.reaction(vcount).reactant(scount).level = 2;
            M.reaction(vcount).reactant(scount).version = 4;
            scount = scount + 1;
            if(~isempty(scomp) && scomp~=ar.model(m).cLink(jsource))
                error('influx from different compartments in reaction %i', jv);
            end
            if(~isempty(ar.model(m).cLink))
                scomp = ar.model(m).cLink(jsource);
            end
        end
        
        %set empty struct for product
        M.reaction(vcount).product = F.reaction(1).product;
        scount = 1;
        tcomp = [];
        for jsource = find(ar.model(m).N(:,jv)>0)'
            M.reaction(vcount).product(scount).typecode = 'SBML_SPECIES_REFERENCE';
            M.reaction(vcount).product(scount).metaid = '';
            M.reaction(vcount).product(scount).notes = '';
            M.reaction(vcount).product(scount).annotation = '';
            M.reaction(vcount).product(scount).sboTerm = -1;
            M.reaction(vcount).product(scount).species = ar.model(m).x{jsource};
            M.reaction(vcount).product(scount).id = '';
            M.reaction(vcount).product(scount).name = '';
            M.reaction(vcount).product(scount).stoichiometry = abs(ar.model(m).N(jsource,jv));
            M.reaction(vcount).product(scount).stoichiometryMath = F.reaction(2).reactant.stoichiometryMath;
            M.reaction(vcount).product(scount).level = 2;
            M.reaction(vcount).product(scount).version = 4;
            scount = scount + 1;
            if(~isempty(tcomp) && tcomp~=ar.model(m).cLink(jsource))
                error('efflux to different compartments in reaction %i', jv);
            end
            if(~isempty(ar.model(m).cLink))
                tcomp = ar.model(m).cLink(jsource);
            end
        end
        
        vars = symvar(ratetemplate);
        vars = setdiff(vars, arSym(ar.model(m).x(ar.model(m).N(:,jv)<0))); %R2013a compatible
        vars = setdiff(vars, arSym(ar.model(m).condition(c).p)); %R2013a compatible
        vars = setdiff(vars, arSym(ar.model(m).u)); %R2013a compatible
        M.reaction(vcount).modifier = F.reaction(1).modifier;
        if(~isempty(vars))
            for jmod = 1:length(vars);
                M.reaction(vcount).modifier(jmod).typecode = 'SBML_MODIFIER_SPECIES_REFERENCE';
                M.reaction(vcount).modifier(jmod).metaid = '';
                M.reaction(vcount).modifier(jmod).notes = '';
                M.reaction(vcount).modifier(jmod).annotation = '';
                M.reaction(vcount).modifier(jmod).sboTerm = -1;
                M.reaction(vcount).modifier(jmod).species = char(vars(jmod));
                M.reaction(vcount).modifier(jmod).id = '';
                M.reaction(vcount).modifier(jmod).name = '';
                M.reaction(vcount).modifier(jmod).level = 2;
                M.reaction(vcount).modifier(jmod).version = 4;
                
            end
        end
        
        M.reaction(vcount).kineticLaw.typecode = 'SBML_KINETIC_LAW';
        M.reaction(vcount).kineticLaw.metaid = '';
        M.reaction(vcount).kineticLaw.notes = '';
        M.reaction(vcount).kineticLaw.annotation = '';
        M.reaction(vcount).kineticLaw.sboTerm = -1;
        
        amountBased = (isfield(ar.model(m),'isAmountBased') && ar.model(m).isAmountBased);
        if ( ~copasi && amountBased )
            warning( 'Exporting amount based models is still in the process of being tested' );
        end
        
        if(~isempty(ar.model(m).cLink))
            if(~copasi)
                M.reaction(vcount).kineticLaw.formula = char(ratetemplate);
            elseif(copasi && amountBased)
                M.reaction(vcount).kineticLaw.formula = [ char(ratetemplate) ]; %' / (' ar.model(m).c{scomp} ')'
            elseif(~isempty(scomp) && ~isempty(tcomp) && scomp~=tcomp) % multi-compartment reaction
                M.reaction(vcount).kineticLaw.formula = [ar.model(m).c{scomp} ' * (' char(ratetemplate) ')'];
            else
                if(~isempty(scomp))
                    M.reaction(vcount).kineticLaw.formula = [ar.model(m).c{scomp} ' * (' char(ratetemplate) ')'];
                elseif(~isempty(tcomp))
                    M.reaction(vcount).kineticLaw.formula = [ar.model(m).c{tcomp} ' * (' char(ratetemplate) ')'];
                else
                    error('scomp and tcomp empty');
                end
            end
        else
            M.reaction(vcount).kineticLaw.formula = char(ratetemplate);
        end
        M.reaction(vcount).kineticLaw.math = M.reaction(vcount).kineticLaw.formula;
        M.reaction(vcount).kineticLaw.parameter = F.reaction(1).kineticLaw.parameter;
        M.reaction(vcount).kineticLaw.level = 2;
        M.reaction(vcount).kineticLaw.version = 4;
        
        vcount = vcount + 1;
    end
end


%% Inputs

% find all possible input functions from arInputFunctionsC.h
% fid = fopen([fileparts(which('arInit')) filesep 'arInputFunctionsC.h'],'r');
% A = fread(fid,'*char')';
% funs = regexp(A,'\ndouble\s(\w*)','tokens')';
% funs = [funs{:}];
funs = cellfun(@(x) x{1}, ar.config.specialFunc,'uniformoutput',0);

for ju = 1:length(ar.model(m).u)
    
    fu = arSym(ar.model(m).condition(c).fu{ju}); % current input
    
    isActive = ~(str2double(ar.model(m).condition(c).fu{ju}) == 0);
    if(contains(ar.model(m).fu{ju},'spline'))
        warning('Spline functions are not supported as d2d export, yet!\n')
        isActive = false; 
    end
    
    if ( isActive )
        % replace p with condition specific parameters
        fu = char(arSubs(fu, arSym(ar.model(m).condition(c).pold), arSym(ar.model(m).condition(c).fp')));
        
        % replace time parameters with 'time'
        fu = char(arSubs(arSym(fu), arSym(ar.model(m).t), arSym('time')));
        ixfun = cell2mat(cellfun(@(x) strfind(fu,x), funs, 'UniformOutput',0)); % does input contain any of the special ar input functions
        if any(ixfun)
            heavisideReplacement = {
                %{'heaviside', 'if((%s) lt 0.0, 0.0, if((%s) gt 0.0, 1.0, 0.5))', [1, 1], 'heaviside(x)'},...
                %{'heaviside', 'if((%s) lt 0.0, 0.0, if((%s) gt 0.0, 1.0, 0.5))', [1, 1], 'heaviside(x)'},...
                {'heaviside', 'piecewise(0, lt((%s), 0), 1)', [1], 'heaviside(x)' }, ...
            };

            % replace functions first
            fu = replaceFunctions( fu, ar.config.specialFunc, 0 );
            % replace heavisides and log their positions
            fu = arSubs(arSym(strrep(fu, 'heaviside', '_potato_')), arSym('_potato_'), 'heaviside');
            [fu, args] = replaceFunctions( fu, heavisideReplacement, 0 );

            % Determine event triggers in here (to make sure SBML resets the solver)
            heaviside_timePoints = cell( numel( args.heaviside ), 1 );
            for jh = 1 : numel( args.heaviside )
                heaviside_timePoints{jh} = args.heaviside(jh).args{1};
            end

            ixrule = length(M.rule) + 1;% index of current rule
            M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
            M.rule(ixrule).metaid = '';
            M.rule(ixrule).notes = '';
            M.rule(ixrule).annotation = '';
            M.rule(ixrule).sboTerm = -1;
            M.rule(ixrule).formula = fu;
            M.rule(ixrule).variable = ar.model(m).u{ju};
            M.rule(ixrule).species = '';
            M.rule(ixrule).compartment = '';
            M.rule(ixrule).name = '';
            if ~isempty(ar.model(m).uUnits)
                M.rule(ixrule).units = ar.model(m).uUnits{ju,2};
            else
                M.rule(ixrule).units = '';
            end
            M.rule(ixrule).level = 2;
            M.rule(ixrule).version = 4;
            if ( 0 )
                for jh = 1 : numel( heaviside_timePoints )
                    %event
                    ixevent = length(M.event) +1;% index of current event
                    M.event(ixevent).typecode =  'SBML_EVENT';
                    M.event(ixevent).metaid = '';
                    M.event(ixevent).notes = '';
                    M.event(ixevent).annotation = '';
                    M.event(ixevent).sboTerm = -1;
                    M.event(ixevent).name = ar.model(m).u{ju};
                    M.event(ixevent).id = sprintf('e%d_%s_event', ixevent, ar.model(m).u{ju});
                    M.event(ixevent).useValuesFromTriggerTime = 1;
                    M.event(ixevent).trigger.typecode =  'SBML_TRIGGER';

                    % construct event trigger
                    % parts = strsplit(fu,{' ',',',')'});
                    M.event(ixevent).trigger.metaid =  '';
                    M.event(ixevent).trigger.notes =  '';
                    M.event(ixevent).trigger.annotation =  '';
                    M.event(ixevent).trigger.sboTerm =  -1;
                    M.event(ixevent).trigger.math = sprintf('gt(%s)',heaviside_timePoints{jh});
                    M.event(ixevent).trigger.level =  2;
                    M.event(ixevent).trigger.version =  4;

                    %         M.event.delay = [1x0 struct];
                    M.event(ixevent).eventAssignment.typecode = 'SBML_EVENT_ASSIGNMENT';
                    M.event(ixevent).eventAssignment.metaid = '';
                    M.event(ixevent).eventAssignment.notes = 'This variable is only used to indicate when events occur so that the solver gets reset appropriately.';
                    M.event(ixevent).eventAssignment.annotation = '';
                    M.event(ixevent).eventAssignment.sboTerm = -1;
                    M.event(ixevent).eventAssignment.variable = 'EventTrigger_dummy'; %ar.model(m).u{ju};
                    M.event(ixevent).eventAssignment.math = '1'; %parts{4};
                    M.event(ixevent).eventAssignment.level = 2;
                    M.event(ixevent).eventAssignment.version = 4;

                    M.event(ixevent).level = 2;
                    M.event(ixevent).version = 4;
                end
            end
        else
            %rule

            ixrule = length(M.rule) +1;% index of current rule

            M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
            M.rule(ixrule).metaid = '';
            M.rule(ixrule).notes = '';
            M.rule(ixrule).annotation = '';
            M.rule(ixrule).sboTerm = -1;
            M.rule(ixrule).formula = fu;
            M.rule(ixrule).variable = ar.model(m).u{ju};
            M.rule(ixrule).species = '';
            M.rule(ixrule).compartment = '';
            M.rule(ixrule).name = '';
            if ~isempty(ar.model(m).uUnits)
                M.rule(ixrule).units = ar.model(m).uUnits{ju,2};
            else
                M.rule(ixrule).units = '';
            end
            M.rule(ixrule).level = 2;
            M.rule(ixrule).version = 4;

        end
    
        initValue = ar.model(m).condition(c).uFineSimu(1,ju);
        % Technically, in the SBML standard, it cannot be considered constant
        % if there is a rule that applies to it. Any value other than zero in
        % this code will fail SBML validation despite being conceptually
        % correct.
        isConstant = 0; %isempty(symvar(fu)); %only cases whith explicit numbers    

        % generate new parameter for each input species
        jp = length(M.parameter)+1;
        M.parameter(jp).typecode = 'SBML_PARAMETER';
        M.parameter(jp).metaid = '';
        M.parameter(jp).notes = '';
        M.parameter(jp).annotation = '';
        M.parameter(jp).sboTerm = -1;
        M.parameter(jp).name = '';
        M.parameter(jp).id = ar.model(m).u{ju};
        if ~isempty(ar.model(m).uUnits)
            M.parameter(jp).units = ar.model(m).uUnits{ju,2};
        else
            M.parameter(jp).units = '';
        end
        M.parameter(jp).constant = isConstant;
        M.parameter(jp).isSetValue = 1;
        M.parameter(jp).level = 2;
        M.parameter(jp).version = 4;
        M.parameter(jp).value = initValue;
    elseif(contains(ar.model(m).fu{ju},'spline_pos'))      
            
            %Getting spline parameters and calculate spline in MATLAB
            nr_ts = strsplit(ar.model.fu{1},'spline_pos');
            nr_ts = strsplit(nr_ts{2},'(');
            nr_ts = str2double(nr_ts{1});
            fprintf('Found cubic interpolation spline with %i anchors! Proceeding with export \n',nr_ts)
            warning('The spline will be re-fitted via a MATLAB function! The result might diverge slightly from the C function implemented in D2D and thus lead to different results when loading the SBML file! \n')
            spline_split = strsplit(ar.model(m).condition(c).fu{ju}, ',');
            spline_times = NaN(1,nr_ts);
            for it = 1:length(spline_times)
                spline_times(it) = str2double(spline_split{it*2});
            end
            par_spline = NaN(1,nr_ts);
            nr_spline = 1;
            for ipu = find(ismember(ar.pLabel,ar.model(m).data(d).pu))
                if ar.qLog10(ipu)
                    par_spline(nr_spline) = 10.^(ar.p(ipu));
                else
                    par_spline(nr_spline) = ar.p(ipu);
                end
                nr_spline = nr_spline + 1;
            end
            par_spline = log(par_spline);           
            pp = csape(spline_times,par_spline,'complete'); %MATLAB solution, nearly as C
            
            if(exist('pp_swameye_Aug18.mat','file')==2)
                load('pp_swameye_Aug18.mat');
                warning('Detected Swameye Export, loaded Spline Parameters');
            end
            %pp.coefs %Coefficient matrix Anchor * nPars  f(x)=a(x?x1)^3+b(x?x1)^2+c(x?x1)+d?.
            %spline_timePoints = pp.breaks; %Anchor points
            
            %Define spline parameters
            for js = 1:(size(pp.coefs,2)+2)
                isConstant = 0;    

                % generate new parameter for each input species
                jp = length(M.parameter)+1;
                M.parameter(jp).typecode = 'SBML_PARAMETER';
                M.parameter(jp).metaid = '';
                M.parameter(jp).notes = '';
                M.parameter(jp).annotation = '';
                M.parameter(jp).sboTerm = -1;
                M.parameter(jp).name = '';
                if(js == 1)
                    M.parameter(jp).id = ar.model(m).u{ju};
                elseif(js == 2)
                    M.parameter(jp).id = 'spline_t_0';
                else
                    M.parameter(jp).id = char(62+js);
                end
                if ~isempty(ar.model(m).uUnits)
                    M.parameter(jp).units = ar.model(m).uUnits{ju,2};
                else
                    M.parameter(jp).units = '';
                end
                M.parameter(jp).constant = isConstant;
                M.parameter(jp).isSetValue = 1;
                
                M.parameter(jp).level = 2;
                M.parameter(jp).version = 4;
                if(js == 1)
                    M.parameter(jp).value = ar.model(m).condition(c).uFineSimu(1,ju);
                else
                    M.parameter(jp).value = 0;
                end
            end
            
            ixrule = length(M.rule) + 1;% index of current rule
            M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
            M.rule(ixrule).metaid = '';
            M.rule(ixrule).notes = '';
            M.rule(ixrule).annotation = '';
            M.rule(ixrule).sboTerm = -1;
            M.rule(ixrule).formula = 'exp(A * (time-spline_t_0)^3 + B * (time-spline_t_0)^2 + C * (time-spline_t_0) + D)';
            M.rule(ixrule).variable = ar.model(m).u{ju};
            M.rule(ixrule).species = '';
            M.rule(ixrule).compartment = '';
            M.rule(ixrule).name = '';
            if ~isempty(ar.model(m).uUnits)
                M.rule(ixrule).units = ar.model(m).uUnits{ju,2};
            else
                M.rule(ixrule).units = '';
            end
            M.rule(ixrule).level = 2;
            M.rule(ixrule).version = 4;
            
            for jh = 1 : (numel(spline_times)-1)
                %event
                ixevent = length(M.event) +1;% index of current event
                M.event(ixevent).typecode =  'SBML_EVENT';
                M.event(ixevent).metaid = '';
                M.event(ixevent).notes = '';
                M.event(ixevent).annotation = '';
                M.event(ixevent).sboTerm = -1;
                M.event(ixevent).name = sprintf('t_%i_event', spline_times(jh));
                M.event(ixevent).id = sprintf('t_%i_event', spline_times(jh));
                M.event(ixevent).useValuesFromTriggerTime = 1;
                M.event(ixevent).level = 2;
                M.event(ixevent).version = 4;
                M.event(ixevent).delay = F.event(1).delay;
                
                % construct event trigger
                % parts = strsplit(fu,{' ',',',')'});
                M.event(ixevent).trigger.typecode =  'SBML_TRIGGER';
                M.event(ixevent).trigger.metaid =  '';
                M.event(ixevent).trigger.notes =  '';
                M.event(ixevent).trigger.annotation =  '';
                M.event(ixevent).trigger.sboTerm =  -1;
                %M.event(ixevent).trigger.math = sprintf('time gt(%s)',num2str(spline_times(jh)));
                M.event(ixevent).trigger.math = sprintf('gt(time,%i)',spline_times(jh));
                M.event(ixevent).trigger.level =  2;
                M.event(ixevent).trigger.version =  4;
                
                %Assign parameters to spline window
                for js = 1:(size(pp.coefs,2)+1)
                    M.event(ixevent).eventAssignment(js).typecode = 'SBML_EVENT_ASSIGNMENT';
                    M.event(ixevent).eventAssignment(js).metaid = '';
                    M.event(ixevent).eventAssignment(js).notes = '';
                    M.event(ixevent).eventAssignment(js).annotation = '';
                    M.event(ixevent).eventAssignment(js).sboTerm = -1;
                    if(js == 1)
                        M.event(ixevent).eventAssignment(js).variable = 'spline_t_0';                   
                        M.event(ixevent).eventAssignment(js).math = mat2str(pp.breaks(jh)); 
                    else
                        M.event(ixevent).eventAssignment(js).variable = char(63+js);                   
                        M.event(ixevent).eventAssignment(js).math = mat2str(pp.coefs(jh,js-1)); 
                    end
                    
                    M.event(ixevent).eventAssignment(js).level = 2;
                    M.event(ixevent).eventAssignment(js).version = 4;
                end
                
            end
        else
            fprintf( 'Input %s not used in condition or is a not supported spline. Omitted from SBML file.\n', ar.model(m).uNames{ju} );            
        end
end
if ( numel(M.event) > 0 )
    jx = numel(M.species) + 1;
    M.species(jx).typecode = 'SBML_SPECIES';
    M.species(jx).metaid = '';
    M.species(jx).notes = '';
    M.species(jx).annotation = '';
    M.species(jx).sboTerm = -1;
    M.species(jx).name = '';
    M.species(jx).id = 'EventTrigger_dummy';
    M.species(jx).speciesType = '';
    M.species(jx).compartment = ar.model(m).c{1};
    M.species(jx).initialAmount = NaN;
    M.species(jx).substanceUnits = '';
    M.species(jx).hasOnlySubstanceUnits = 0;
    M.species(jx).boundaryCondition = 0;
    M.species(jx).charge = 0;
    M.species(jx).constant = 0;
    M.species(jx).isSetInitialAmount = 0;
    M.species(jx).isSetInitialConcentration = 1;
    M.species(jx).isSetCharge = 0;
    M.species(jx).level = 2;
    M.species(jx).version = 4;
    M.species(jx).initialConcentration = 0;
end

% assign time symbol
%M.time_symbol = ar.model(m).t;

M.unitDefinition = F.unitDefinition;
if(strcmp(ar.model(m).tUnits{2},'h'))
    M.unitDefinition(1).unit(1).multiplier = 3600;
elseif(strcmp(ar.model(m).tUnits{2},'s') || strcmp(ar.model(m).tUnits{2},'sec'))
    M.unitDefinition(1).unit(1).multiplier = 1;
end


%% Observables as assignment rules

for id = 1:length(ar.model(m).data(d).fy)
    ixrule = length(M.rule) + 1;% index of current rule
    M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
    M.rule(ixrule).metaid = '';
    M.rule(ixrule).notes = '';
    M.rule(ixrule).annotation = '';
    M.rule(ixrule).sboTerm = -1;
    rule_tmp = arSym(ar.model(m).data(d).fy{id});
    if(~isempty(ar.model(m).z))
        rule_tmp = arSubs(rule_tmp,arSym(ar.model(m).z),arSym(ar.model(m).fz'));
    end
    rule_tmp = arSubs(rule_tmp,arSym(ar.model(m).data(d).pold),arSym(ar.model(m).data(d).fp'));
    
    rule_tmp = char(rule_tmp);
    if(ar.model(m).data(d).logfitting(id)==1)
        rule_tmp = ['log10(' rule_tmp ')'];
    end
    
    % check if ys already contain 'observable_'
    if ~strncmp(ar.model(m).data(d).y{id}, 'observable_', length('observable_'))
        variable_tmp = strcat('observable_', ar.model(m).data(d).y{id});
    else
        variable_tmp = ar.model(m).data(d).y{id};
    end
    
    M.rule(ixrule).formula = rule_tmp;
    M.rule(ixrule).variable = variable_tmp;
    M.rule(ixrule).species = '';
    M.rule(ixrule).compartment = '';
    M.rule(ixrule).name = '';
    if ~isempty(ar.model(m).data(d).yUnits)
        M.rule(ixrule).units = ar.model(m).data(d).yUnits{id,2};
    else
        M.rule(ixrule).units = '';
    end
    M.rule(ixrule).level = 2;
    M.rule(ixrule).version = 4;
end

%% Errors as assignment rules

for id = 1:length(ar.model(m).data(d).fystd)
    ixrule = length(M.rule) + 1;% index of current rule
    M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
    M.rule(ixrule).metaid = '';
    M.rule(ixrule).notes = '';
    M.rule(ixrule).annotation = '';
    M.rule(ixrule).sboTerm = -1;
    rule_tmp = arSym(ar.model(m).data(d).fystd{id});
    
    rule_tmp = arSubs(rule_tmp,arSym(ar.model(m).data(d).pold),arSym(ar.model(m).data(d).fp'));
    
    rule_tmp = char(rule_tmp);
    
    % ??
    if(ar.model(m).data(d).logfitting(id)==1)
        rule_tmp = ['log10(' rule_tmp ')'];
    end
    
    % check if ys already contain 'observable_'
    if ~strncmp(ar.model(m).data(d).y{id}, 'observable_', length('observable_'))
        variable_tmp = strcat('sigma_', ar.model(m).data(d).y{id});
    else
        variable_tmp = strrep(ar.model(m).data(d).y{id}, 'observable_', 'sigma_');
    end
    M.rule(ixrule).variable = variable_tmp;
    M.rule(ixrule).formula = rule_tmp;
    M.rule(ixrule).species = '';
    M.rule(ixrule).compartment = '';
    M.rule(ixrule).name = '';
    if ~isempty(ar.model(m).data(d).yUnits)
        M.rule(ixrule).units = ar.model(m).data(d).yUnits{id,2};
    else
        M.rule(ixrule).units = '';
    end
    M.rule(ixrule).level = 2;
    M.rule(ixrule).version = 4;
end

arWaitbar(-1);

[a,b] = isSBML_Model(M);
if(a == 1)
%     if(~copasi)
%         OutputSBML(M, ['SBML/' ar.model(m).name '_cond_' num2str(c) '_l2v4.xml']);
%     else
%         OutputSBML(M, ['SBML/' ar.model(m).name '_cond_' num2str(c)  '_copasi_l2v4.xml']);
%     end
    
%mine!
                    if(~exist('./Benchmark_paper/SBML', 'dir'))
                    mkdir('./Benchmark_paper/SBML')
                    end
                    OutputSBML(M, ['Benchmark_paper/SBML/model' num2str(m) '_data' num2str(d) '_l2v4.xml']);
else
    error('%s', b);
end





