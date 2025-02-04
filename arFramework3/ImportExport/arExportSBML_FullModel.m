% function arExportSBML_FullModel(m)
%
% formerly arExportSBML_benchmark
%
% export single data file of model to SBML using libSBML
%
% m:            model index

function arExportSBML_FullModel(m,name)

    global ar
    
    % simulate once for initial values
    try
        arSimu(0,1,0)
    catch
        arSimu(1,1,0)
    end
    
    copasi = true;
    
    % if(~exist('data','var'))
    %     data = 1;
    % end
    
    % if (~exist('steadystate','var'))
    %     if ( isfield( ar.model(m), 'ss_condition' ) )
    %         steadystate = true;
    %     else
    %         steadystate = false;
    %     end
    % end

    qCondsSBMLConform = arCheckSBMLCompatibility(m);
    if ~qCondsSBMLConform
        warning(['Model conditions are not independet and will be represented incorrectly in SBML. ', ...
                 'Consider using "arRenameModelCondPars" to get independent model parameters and conditions.']);
    end
    
    try
        M = TranslateSBML(which('empty.xml'));
        F = TranslateSBML(which('filled.xml'));
    catch
        warning('error in libSBML. Probably backwards compatibility issues with old MATLAB version. Should work with 2019a.')
    end
    
    M.id = ar.model(m).name;
    if(~isempty(ar.model(m).description))
        M.notes = ar.model(m).description{1};
    end
    
    
    [M] = GetCompartments(M,m);
    
    [M,Crules] = GetSpecies(M,m);
    
    [M] = GetParameters(M,m);
    
    [M] = GetInitialAssignments(M,Crules);
    
    [M] = GetReactions(M,F,m,copasi);
    
    [M] = GetInputs(M,F,m);
    
    [M] = GetUnitDefinitions(M,m);
    
    
    % [M] = GetObservables(M,m);
    %
    % [M] = GetErrors(M,m);
    
    
    arWaitbar(-1);
    
    [a,b] = isSBML_Model(M);
    if(a == 1)
        
        %     if(~copasi)
        %         OutputSBML(M, ['SBML/' ar.model(m).name '_cond_' num2str(c) '_l2v4.xml']);
        %     else
        %         OutputSBML(M, ['SBML/' ar.model(m).name '_cond_' num2str(c)  '_copasi_l2v4.xml']);
        %     end
        % mine!
        
        if(~exist('./PEtab', 'dir'))
            mkdir('./PEtab')
        end
        if length(ar.model)==1
            OutputSBML(M, ['PEtab/' name '_model.xml']);
        else
            OutputSBML(M, ['PEtab/' name '_' ar.model(m).name '_model.xml']);
        end
    else
        error('%s', b);
    end
    
    
    end
    
    
    
    %% compartements
    function [M] = GetCompartments(M,m)
        global ar
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
                M.compartment(jc).isSetSize = 0;
                M.compartment(jc).isSetVolume = 0;
                M.compartment(jc).level = 2;
                M.compartment(jc).version = 4;
                M.compartment(jc).size = 1;
                
                % Get the size of the compartment
                cSize = ar.model(m).pc{jc};
                cSizeNum = str2num(cSize);

                if ~isnan(cSizeNum)
                    % Case 1: compartment size is a numeric value
                    % -> directly set the size in the compartment definition
                    M.compartment(jc).size = cSizeNum;
                    M.compartment(jc).isSetSize = 1;
                    M.compartment(jc).isSetVolume = 1;
                else
                    % Cases 2x: compartment size is a parameter (+ a replacement by conditions)
                    % -> add an initial assignment rule to sbml file

                    if ismember(ar.pLabel, cSize)
                        % Case 2a: compartment size is an optimization parameter
                        % -> replace parameter by numerical value from ar.p
                        idx = find(strcmp(ar.pLabel, cSize));
                        cSizeNum = ar.p(idx);
                        if(ar.qLog10(idx))
                            cSizeNum = 10^cSizeNum;
                        end
                        cSizeReplace = num2str(cSizeNum);

                    else
                        % Case 2b/c: compartment size is not an optimization parameter
                        % -> it is only a model parameter (with replacements ind model.def or data.def files)
                        idx = find(strcmp(ar.model(m).p, cSize));
                        modelfp = string(arSym(ar.model(m).fp{idx}));
                        if ~strcmp(modelfp, cSize)
                            % Case 2b: compartment size is replaced in model condition
                            % -> replace parameter accordingly
                            cSizeReplace = modelfp;
                            % cSizeReplace = cSize;

                        else
                            % Case 2c: no model-wide replacement found
                            % -> add the parameter itself to the sbml file
                            cSizeReplace = cSize;

                            % consistency check: Is compartment size replaced in all data conditions?
                            qDataRepl = false(1, length(ar.model(m).data));
                            for d = 1:length(ar.model(m).data)
                                idx = find(strcmp(ar.model(m).data(d).pold, cSize));
                                datafp = string(arSym(ar.model(m).data(d).fp{idx}));
                                qDataRepl(d) = ~strcmp(datafp, cSize);
                            end
                            if ~all(qDataRepl)
                                warning(['Compartment size %s is not defined properly.\n' ...
                                    'Add a numeric value, a condition in %s.def or conditions in all Data/*.def files.'], ...
                                    cSize, ar.model(m).name);
                            end
                            
                        end
                    end

                    % Add an initial assignment rule for the compartment size
                    idxIA = length(M.initialAssignment) + 1;
                    M.initialAssignment(idxIA).typecode = 'SBML_INITIAL_ASSIGNMENT';
                    M.initialAssignment(idxIA).metaid = '';
                    M.initialAssignment(idxIA).notes = '';
                    M.initialAssignment(idxIA).annotation = '';
                    M.initialAssignment(idxIA).sboTerm = -1;
                    M.initialAssignment(idxIA).symbol = M.compartment(jc).id;
                    M.initialAssignment(idxIA).math = char(cSizeReplace);
                    M.initialAssignment(idxIA).level = 2;
                    M.initialAssignment(idxIA).version = 4;

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
    end
    
    
    
    function [M,Crules] = GetSpecies(M,m)
    global ar
    
    %% species
    Crules = {};
    
    
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
               
        simulated_ss = 0;
        %         if ( steadystate )
        %             if (isfield(ar.model(m).condition(c),'ssLink') && ~isempty( ar.model(m).condition(c).ssLink ) )
        %                 warning( 'Using simulated steady state values as initial condition (non-parametric)' );
        %                 x_ss = ar.model(m).ss_condition(ar.model(m).condition(c).ssLink).xFineSimu(end,:);
        %                 simulated_ss = 1;
        %             end
        %         end
        
        if ( ~simulated_ss )
            % check if init parameter still exists in condition parameters
            qp = ismember(ar.pLabel, ar.model(m).px0{jx}); %R2013a compatible
            is_set = sum(ismember(ar.model(m).fp, ar.model(m).px0{jx}))==0;
            if(sum(qp)==1 && is_set==0)
                M.species(jx).initialConcentration = 1;
                Crules{end+1,1} = ar.model(m).x{jx}; %#ok<AGROW>
                Crules{end,2} = ar.pLabel{qp}; %#ok<AGROW>
            elseif(sum(qp)==0 || is_set)
                qp = ismember(ar.model(m).p, ar.model(m).px0{jx}); %R2013a compatible
                if(sum(qp)==1)
                    pvalue = char(arSym(ar.model(m).fp{qp}));
                    % pvalue = char(arSym(ar.model(m).p{qp}));
                    %                     if(~isnan(str2num(pvalue))) %#ok
                    %                         pvalue = str2num(pvalue); %#ok
                    %                         M.species(jx).initialConcentration = pvalue;
                    %                     else
                    Crules{end+1,1} = ar.model(m).x{jx}; %#ok<AGROW>
                    Crules{end,2} = pvalue; %#ok<AGROW>
                    M.species(jx).initialConcentration = 1;
                    %                     end
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
    end
    
    
    
    function [M] = GetUnitDefinitions(M,m)
    global ar
    
    % collect all units
    allUnits = {};
    fields = {ar.model(m).xUnits, ar.model(m).uUnits, ar.model(m).cUnits};
    if ~isempty(ar.model(m).yUnits)
        fields{end+1} = ar.model(m).yUnits;
    end
    for id = 1:length(ar.model(m).data)
        fields{end+1} = ar.model(m).data(id).yUnits;
    end
    for jf = 1:length(fields)
        if ~isempty(fields{jf})
            for jfi = 1:size(fields{jf},1)
                allUnits{end+1} = fields{jf}{jfi,2};
            end
        end
    end
    allUnits = unique(allUnits);
    
    for jun = 1:length(allUnits)
        theUnit = allUnits{jun};
        
        if strcmp(theUnit, 'pl')
            myname = 'pl';
            mykind = 'litre';
            mymultiplier = 1e-12;
        elseif strcmp(theUnit, 'nM')
            myname = 'nM';
            mykind = 'mole';
            mymultiplier = 1e-9;
        else
            warning(sprintf('Unknown unit %s, may be incorrectly parsed to SBML', theUnit))
            continue
        end
        
        ixrule = length(M.unitDefinition) + 1;% index of current rule
        
        M.unitDefinition(ixrule).typecode = 'SBML_UNIT_DEFINITION';
        M.unitDefinition(ixrule).metaid = '';
        M.unitDefinition(ixrule).notes = '';
        M.unitDefinition(ixrule).annotation = '';
        M.unitDefinition(ixrule).cvterms = [];
        M.unitDefinition(ixrule).sboTerm = -1;
        M.unitDefinition(ixrule).name = myname;
        M.unitDefinition(ixrule).id = myname;
        
        M.unitDefinition(ixrule).unit.typecode = 'SBML_UNIT';
        M.unitDefinition(ixrule).unit.metaid = '';
        M.unitDefinition(ixrule).unit.notes = '';
        M.unitDefinition(ixrule).unit.annotation = '';
        M.unitDefinition(ixrule).unit.cvterms = [];
        M.unitDefinition(ixrule).unit.sboTerm = -1;
        M.unitDefinition(ixrule).unit.kind = mykind;
        M.unitDefinition(ixrule).unit.exponent = 1;
        M.unitDefinition(ixrule).unit.scale = 0;
        M.unitDefinition(ixrule).unit.multiplier = mymultiplier;
        M.unitDefinition(ixrule).unit.level = 2;
        M.unitDefinition(ixrule).unit.version = 4;
        
        M.unitDefinition(ixrule).level = 2;
        M.unitDefinition(ixrule).version = 4;
    end
    end
    
    
    function [M] = GetParameters(M,m)
    global ar


    %% first: collect available numerical values for model parameters (ar.model.p)

    % logical flags for model parameters
    isInit = cellfun(@(x) any(strcmp(x, ar.model(m).px0)), ar.model(m).p)';
    isReplaced = ~strcmp(ar.model.p', string(arSym(ar.model.fp)));
    isCompSize = cellfun(@(x) any(strcmp(x, string(arSym(ar.model(m).pc)))), ar.model(m).p)';
    % isConst = (isReplaced & ~isInit) | (isInit & ~isReplaced);
    isConst = true(1, length(ar.model(m).p));  % all model parameters are constant (i.e. not time-dependent)

    for jp = 1:length(ar.model(m).p)
        % All model parameters should be defined in SBML as parameters
        % irrespective of initAssigns or undefined values.

        % compound parameters are handeled separately in "GetCompartments"
        if isCompSize(jp)
            continue
        end

        % Is there a numeric value for the model parameter in ar.p?
        qp = strcmp(ar.model(m).p(jp), ar.pLabel); %R2013a compatible
        if any(qp)
            % get parameter value from ar.p
            pvalue = ar.p(qp);
            if(ar.qLog10(qp) == 1)
                pvalue = 10^pvalue;
            end
            isSetValue = 1;
            constant = double(isConst(qp));
        else
            % no numeric value in ar.p found
            pvalue = NaN;
            isSetValue = 0;
            constant = 1;
        end

        id_tmp = length(M.parameter) + 1;
        M.parameter(id_tmp).typecode = 'SBML_PARAMETER';
        M.parameter(id_tmp).metaid = '';
        M.parameter(id_tmp).notes = '';
        M.parameter(id_tmp).annotation = '';
        M.parameter(id_tmp).sboTerm = -1;
        M.parameter(id_tmp).name = ar.model(m).p{jp};
        M.parameter(id_tmp).id = ar.model(m).p{jp};
        M.parameter(id_tmp).units = '';
        M.parameter(id_tmp).constant = constant;
        M.parameter(id_tmp).isSetValue = isSetValue;
        M.parameter(id_tmp).value = pvalue;
        M.parameter(id_tmp).level = 2;
        M.parameter(id_tmp).version = 4;
    
    end


    %% second: collect replacements of model parameters

    for jp = 1:length(ar.model(m).p)
        if isReplaced(jp) && ~isCompSize(jp)

            if isInit(jp)
                % initial values for species are already defined in "GetSpecies"
                % and written to SBML in "GetInitialAssignments".
                % However: An init parameter can also appears as a model parameter
                % (e.g. in a rate equation or input). Then the CONDITION
                %   init_State  "expression"
                % must also be applied to the parameter, not just the species.
                if ~ismember(ar.model(m).p{jp}, union(ar.model(m).pu, ar.model(m).pv))
                    continue
                end
            end
            
            % parameter CONDITIONS should be implemeted as initialAssignment
            % reason: parameters in d2d are constant, CONDITIONS are applied before start of simulation

            ixInitAssign = length(M.initialAssignment) + 1;% index of current rule
            M.initialAssignment(ixInitAssign).typecode = 'SBML_INITIAL_ASSIGNMENT';
            M.initialAssignment(ixInitAssign).metaid = '';
            M.initialAssignment(ixInitAssign).notes = '';
            M.initialAssignment(ixInitAssign).annotation = '';
            M.initialAssignment(ixInitAssign).sboTerm = -1;
            M.initialAssignment(ixInitAssign).symbol = ar.model(m).p{jp};
            M.initialAssignment(ixInitAssign).math = ar.model(m).fp{jp};
            M.initialAssignment(ixInitAssign).level = 2;
            M.initialAssignment(ixInitAssign).version = 4;

            % ixrule = length(M.rule) + 1;% index of current rule
            % M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
            % M.rule(ixrule).metaid = '';
            % M.rule(ixrule).notes = '';
            % M.rule(ixrule).annotation = '';
            % M.rule(ixrule).sboTerm = -1;
            % M.rule(ixrule).formula = ar.model(m).fp{jp};
            % M.rule(ixrule).variable = ar.model(m).p{jp};
            % M.rule(ixrule).species = '';
            % M.rule(ixrule).compartment = '';
            % M.rule(ixrule).name = '';
            % M.rule(ixrule).units = '';
            % M.rule(ixrule).level = 2;
            % M.rule(ixrule).version = 4;
        end    
    end
    
    %% third: search numerical values for replacements
    
    % optimization parameters that appear in replacements (and are not model parameters)
    replParams = cellfun(@(x) any(contains(ar.model(m).fp(isReplaced), x)), ar.pLabel);
    replParams = replParams & ~cellfun(@(x) ismember(x, ar.model(m).p), ar.pLabel);
    
    for jp = 1:length(ar.pLabel)
        
        if replParams(jp)
   
            id_tmp = length(M.parameter) + 1;
            M.parameter(id_tmp).typecode = 'SBML_PARAMETER';
            M.parameter(id_tmp).metaid = '';
            M.parameter(id_tmp).notes = '';
            M.parameter(id_tmp).annotation = '';
            M.parameter(id_tmp).sboTerm = -1;
            M.parameter(id_tmp).name = ar.pLabel{jp};
            M.parameter(id_tmp).id = ar.pLabel{jp};
            M.parameter(id_tmp).units = '';
            M.parameter(id_tmp).constant = 1;
            M.parameter(id_tmp).isSetValue = 1;
            M.parameter(id_tmp).level = 2;
            M.parameter(id_tmp).version = 4;
            
            pvalue = ar.p(jp);
            if(ar.qLog10(jp) == 1)
                pvalue = 10^pvalue;
            end
            M.parameter(id_tmp).value = pvalue;
            
        end
    end
    end
    
    
    function [M] = GetInitialAssignments(M,Crules)
    
    %% rules (copasi)
    idxStartIA = length(M.initialAssignment);
    for jr = 1:size(Crules,1)
        idxIA = idxStartIA + jr;
        M.initialAssignment(idxIA).typecode = 'SBML_INITIAL_ASSIGNMENT';
        M.initialAssignment(idxIA).metaid = '';
        M.initialAssignment(idxIA).notes = '';
        M.initialAssignment(idxIA).annotation = '';
        M.initialAssignment(idxIA).sboTerm = -1;
        M.initialAssignment(idxIA).symbol = Crules{jr,1};
        M.initialAssignment(idxIA).math = Crules{jr,2};
        M.initialAssignment(idxIA).level = 2;
        M.initialAssignment(idxIA).version = 4;
    end
    
    end
    
    
    
    
    function [M] = GetReactions(M,F,m,copasi)
    global ar
    %% reactions
    fv = ar.model(m).fv;
    fv = arSym(fv);
    % fv = arSubs(fv, ar.model(m).u, ar.model(m).condition(c).fu');
    %fv = arSubs(fv, arSym(ar.model(m).condition(c).pold), arSym(ar.model(m).condition(c).fp'));
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
            vars = setdiff(vars, arSym(ar.model(m).p)); %R2013a compatible
            vars = setdiff(vars, arSym(ar.model(m).u)); %R2013a compatible
            M.reaction(vcount).modifier = F.reaction(1).modifier;
            if(~isempty(vars))
                for jmod = 1:length(vars)
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
    
    
    end
    
    
    
    
    function [M] = GetInputs(M,F,m)
    global ar
    %% Inputs
    
    % find all possible input functions from arInputFunctionsC.h
    % fid = fopen([fileparts(which('arInit')) filesep 'arInputFunctionsC.h'],'r');
    % A = fread(fid,'*char')';
    % funs = regexp(A,'\ndouble\s(\w*)','tokens')';
    % funs = [funs{:}];
    funs = cellfun(@(x) x{1}, ar.config.specialFunc,'uniformoutput',0);
    
    for ju = 1:length(ar.model(m).u)
        
        fu = arSym(ar.model(m).fu{ju});
        
        isActive = ~(str2double(ar.model(m).fu{ju}) == 0);
        if(contains(ar.model(m).fu{ju},'spline'))
            warning('Spline functions are not supported as d2d export, yet!\n')
            isActive = false;
        end
        
        if ( isActive )
            
            % replace time parameters with 'time'
            fu = char(arSubs(arSym(fu), arSym(ar.model(m).t), arSym('time')));
            ixfun = cell2mat(cellfun(@(x) strfind(fu,x), funs, 'UniformOutput',0)); % does input contain any of the special ar input functions
            if any(ixfun)
                % because subs will be called time has to be renamed
                fu = arSym(ar.model(m).fu{ju});
                fu = char(arSubs(arSym(fu), arSym(ar.model(m).t), arSym('cabbage_')));
                
                heavisideReplacement = {
                    %{'heaviside', 'if((%s) lt 0.0, 0.0, if((%s) gt 0.0, 1.0, 0.5))', [1, 1], 'heaviside(x)'},...
                    %{'heaviside', 'if((%s) lt 0.0, 0.0, if((%s) gt 0.0, 1.0, 0.5))', [1, 1], 'heaviside(x)'},...
                    {'heaviside', 'piecewise(0, lt((%s), 0), 1)', [1], 'heaviside(x)' }, ...
                    };
                
                % replace functions first
                fu = replaceFunctions( fu, ar.config.specialFunc, 0 );
                % replace heavisides and log their positions
                fu = arSubs(arSym(strrep(fu, 'heaviside', 'potato_')), arSym('potato_'), arSym('heaviside'));
                [fu, args] = replaceFunctions( fu, heavisideReplacement, 0 );
                
                % Determine event triggers in here (to make sure SBML resets the solver)
                heaviside_timePoints = cell( numel( args.heaviside ), 1 );
                for jh = 1 : numel( args.heaviside )
                    heaviside_timePoints{jh} = args.heaviside(jh).args{1};
                end
                
                fu = strrep(fu, 'cabbage_', 'time');
                
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
                        M.event(ixevent).trigger.math = sprintf('ge(%s)',heaviside_timePoints{jh});
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
            
            %initValue = ar.model(m).condition(c).uFineSimu(1,ju);
            initValue = 0;
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
            M.parameter(jp).name = ar.model(m).u{ju};
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
            spline_split = strsplit(ar.model(m).fu{ju}, ',');
            spline_times = NaN(1,nr_ts);
            for it = 1:length(spline_times)
                spline_times(it) = str2double(spline_split{it*2});
            end
            par_spline = NaN(1,nr_ts);
            nr_spline = 1;
            for ipu = find(ismember(ar.pLabel,ar.model(m).pu))
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
                
                M.parameter(jp).value = 0;
                
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
                M.event(ixevent).trigger.math = sprintf('ge(time,%i)',spline_times(jh));
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
            fprintf( 'Input %s not used in condition or is a not supported spline. Omitted from SBML file.\n', ar.model(m).u{ju} );
            if ~isnan(str2double(ar.model(m).fu{ju}))
                warning('Input %s is a numeric value. Please use step1-function instead of integer in model definition to enable export to SBML format.\n', ar.model(m).u{ju} );
            end
        end
    end
    
    %% Species and units
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
    
    end
    
    
    
    
    function [M] = GetObservables(M,m)
    global ar
    %% Observables as assignment rules
    
    for id = 1:length(ar.model(m).fy)
        ixrule = length(M.rule) + 1;% index of current rule
        M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
        M.rule(ixrule).metaid = '';
        M.rule(ixrule).notes = '';
        M.rule(ixrule).annotation = '';
        M.rule(ixrule).sboTerm = -1;
        rule_tmp = arSym(ar.model(m).fy{id});
        if(~isempty(ar.model(m).z))
            rule_tmp = arSubs(rule_tmp,arSym(ar.model(m).z),arSym(ar.model(m).fz'));
        end
        %rule_tmp = arSubs(rule_tmp,arSym(ar.model(m).data(d).pold),arSym(ar.model(m).data(d).fp'));
        
        rule_tmp = char(rule_tmp);
        % logfitting is specfied in peTAB obs file!
        %         if(ar.model(m).logfitting(id)==1)
        %             rule_tmp = ['log10(' rule_tmp ')'];
        %         end
        
        % check if ys already contain 'observable_'
        if ~strncmp(ar.model(m).y{id}, 'observable_', length('observable_'))
            variable_tmp = strcat('observable_', ar.model(m).y{id});
        else
            variable_tmp = ar.model(m).y{id};
        end
        
        M.rule(ixrule).formula = rule_tmp;
        M.rule(ixrule).variable = variable_tmp;
        M.rule(ixrule).species = '';
        M.rule(ixrule).compartment = '';
        M.rule(ixrule).name = '';
        if ~isempty(ar.model(m).yUnits)
            M.rule(ixrule).units = ar.model(m).yUnits{id,2};
        else
            M.rule(ixrule).units = '';
        end
        M.rule(ixrule).level = 2;
        M.rule(ixrule).version = 4;
    end
    
    
    
    end
    
    
    
    function [M] = GetErrors(M,m)
    global ar
    %% Errors as assignment rules
    
    for id = 1:length(ar.model(m).fystd)
        ixrule = length(M.rule) + 1;% index of current rule
        M.rule(ixrule).typecode = 'SBML_ASSIGNMENT_RULE';
        M.rule(ixrule).metaid = '';
        M.rule(ixrule).notes = '';
        M.rule(ixrule).annotation = '';
        M.rule(ixrule).sboTerm = -1;
        rule_tmp = arSym(ar.model(m).fystd{id});
        
        %rule_tmp = arSubs(rule_tmp,arSym(ar.model(m).data(d).pold),arSym(ar.model(m).data(d).fp'));
        
        rule_tmp = char(rule_tmp);
        
        % logfitting is specfied in peTAB obs file!
        %         if(ar.model(m).logfitting(id)==1)
        %             rule_tmp = ['log10(' rule_tmp ')'];
        %         end
        
        % check if ys already contain 'observable_'
        if ~strncmp(ar.model(m).y{id}, 'observable_', length('observable_'))
            variable_tmp = strcat('sigma_', ar.model(m).y{id});
        else
            variable_tmp = strrep(ar.model(m).y{id}, 'observable_', 'sigma_');
        end
        M.rule(ixrule).variable = variable_tmp;
        M.rule(ixrule).formula = rule_tmp;
        M.rule(ixrule).species = '';
        M.rule(ixrule).compartment = '';
        M.rule(ixrule).name = '';
        if ~isempty(ar.model(m).yUnits)
            M.rule(ixrule).units = ar.model(m).yUnits{id,2};
        else
            M.rule(ixrule).units = '';
        end
        M.rule(ixrule).level = 2;
        M.rule(ixrule).version = 4;
    end
    
    
    end
    
    