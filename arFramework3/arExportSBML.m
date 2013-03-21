% export model to SBML using libSBML
%
% function arExportSBML(m, c, copasi)
%
% m:        model index
% c:        condition index
% copasi:   copasi compatibility mode

function arExportSBML(m, c, copasi)

global ar

if(~exist([cd '/SBML' ], 'dir'))
    mkdir([cd '/SBML' ])
end

if(~exist('copasi','var'))
    copasi = false;
end

M = TranslateSBML(which('empty.xml'));
F = TranslateSBML(which('filled.xml'));

M.id = ar.model(m).name;
M.notes = ar.model(m).description{1};

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
        M.compartment(jc).units = '';
        M.compartment(jc).outside = '';
        M.compartment(jc).isSetSize = 1;
        M.compartment(jc).isSetVolume = 1;
        M.compartment(jc).level = 2;
        M.compartment(jc).version = 4;
        
        qp = ismember(ar.pLabel, ar.model(m).pc{jc});
        if(sum(qp)==1)
            pvalue = ar.p(qp);
            if(ar.qLog10(qp))
                pvalue = 10^pvalue;
            end
            M.compartment(jc).size = pvalue;
        elseif(sum(qp)==0)
            qp = ismember(ar.model(m).condition(c).pold, ar.model(m).pc{jc});
            if(sum(qp)==1)
                pvalue = ar.model(m).condition(c).fp{qp};
                M.compartment(jc).size = str2double(pvalue);
            else
                error('%s not found', ar.model(m).pc{jc});
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
    M.compartment(1).units = '';
    M.compartment(1).outside = '';
    M.compartment(1).isSetSize = 1;
    M.compartment(1).isSetVolume = 1;
    M.compartment(1).level = 2;
    M.compartment(1).version = 4;
    M.compartment(1).size = 1;
end

%% species
Crules = {};
for jx = 1:length(ar.model(m).x)
    
    M.species(jx).typecode = 'SBML_SPECIES';
    M.species(jx).metaid = '';
    M.species(jx).notes = '';
    M.species(jx).annotation = '';
    M.species(jx).sboTerm = -1;
    M.species(jx).name = '';
    M.species(jx).id = ar.model(m).x{jx};
    M.species(jx).speciesType = '';
    if(~isempty(ar.model(m).cLink))
        M.species(jx).compartment = ar.model(m).c{ar.model(m).cLink(jx)};
    else
        M.species(jx).compartment = 'default';
    end
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
    
    qp = ismember(ar.pLabel, ar.model(m).px0{jx});
    if(sum(qp)==1)
        M.species(jx).initialConcentration = 1;
        Crules{end+1,1} = ar.model(m).x{jx}; %#ok<AGROW>
        Crules{end,2} = ar.pLabel{qp}; %#ok<AGROW>
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).px0{jx});
        if(sum(qp)==1)
            pvalue = char(sym(ar.model(m).condition(c).fp{qp}));
            if(~isnan(str2double(pvalue)))
                pvalue = str2double(pvalue);
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
end

%% parameters
for jp = 1:length(ar.model(m).condition(c).p)
    M.parameter(jp).typecode = 'SBML_PARAMETER';
    M.parameter(jp).metaid = '';
    M.parameter(jp).notes = '';
    M.parameter(jp).annotation = '';
    M.parameter(jp).sboTerm = -1;
    M.parameter(jp).name = '';
    M.parameter(jp).id = ar.model(m).condition(c).p{jp};
    M.parameter(jp).units = '';
    M.parameter(jp).constant = 1;
    M.parameter(jp).isSetValue = 1;
    M.parameter(jp).level = 2;
    M.parameter(jp).version = 4;
    
    qp = ismember(ar.pLabel, ar.model(m).condition(c).p{jp});
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp) == 1)
            pvalue = 10^pvalue;
        end
    else
        pvalue = 1;
    end
    M.parameter(jp).value = pvalue;
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
fv = sym(fv);
fv = subs(fv, ar.model(m).u, ar.model(m).fu);
fv = subs(fv, ar.model(m).condition(c).pold, ar.model(m).condition(c).fp);

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
        M.reaction(vcount).name = '';
        M.reaction(vcount).reversible = 0;
        M.reaction(vcount).fast = -1;
        M.reaction(vcount).isSetFast = 0;
        M.reaction(vcount).level = 2;
        M.reaction(vcount).version = 4;
        
        if(isfield(ar.model(m),'v') && ~isempty(ar.model(m).v{jv}))
            M.reaction(vcount).id = ar.model(m).v{jv};
        else
            M.reaction(vcount).id = sprintf('reaction%i', jv);
        end
        
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
            M.reaction(vcount).reactant(scount).stoichiometryMath = F.reaction.reactant.stoichiometryMath;
            M.reaction(vcount).reactant(scount).level = 2;
            M.reaction(vcount).reactant(scount).version = 4;
            scount = scount + 1;
            if(~isempty(scomp) && scomp~=ar.model.cLink(jsource))
                error('influx from different compartments in reaction %i', jv);
            end
            if(~isempty(ar.model(m).cLink))
                scomp = ar.model.cLink(jsource);
            end
        end
        
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
            M.reaction(vcount).product(scount).stoichiometryMath = F.reaction.reactant.stoichiometryMath;
            M.reaction(vcount).product(scount).level = 2;
            M.reaction(vcount).product(scount).version = 4;
            scount = scount + 1;
            if(~isempty(tcomp) && tcomp~=ar.model.cLink(jsource))
                error('efflux to different compartments in reaction %i', jv);
            end
            if(~isempty(ar.model(m).cLink))
                tcomp = ar.model.cLink(jsource);
            end
        end
        
        vars = symvar(ratetemplate);
        vars = setdiff(vars, sym(ar.model(m).x(ar.model(m).N(:,jv)<0)));
        vars = setdiff(vars, sym(ar.model(m).condition(c).p));
        
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
        if(~isempty(ar.model(m).cLink))
            if(~copasi)
                M.reaction(vcount).kineticLaw.formula = char(ratetemplate);
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
        M.reaction(vcount).kineticLaw.math = '';
        M.reaction(vcount).kineticLaw.parameter = F.reaction.kineticLaw.parameter;
        M.reaction(vcount).kineticLaw.level = 2;
        M.reaction(vcount).kineticLaw.version = 4;
        
        vcount = vcount + 1;
    end
end
arWaitbar(-1);

[a,b] = isSBML_Model(M);
if(a == 1)
    OutputSBML(M, ['SBML/' ar.model(m).name '_l2v4.xml']);
else
    error('%s', b);
end


