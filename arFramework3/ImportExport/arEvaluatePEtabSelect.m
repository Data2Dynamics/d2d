% status = arEvaluatePEtabSelect(modelSpecFile, constrFile, config)
% 
% arEvaluatePEtabSelect calculates scores like AIC, BIC, etc for a set of
% candidate models. 
% 
%   modelSpecFile   *.tsv file that contains info about all candidate models that
%                   should be calculated
%   [constrFile]    Additional constraints that should be fulfilled when
%                   investigating candidate models
% 
%   status          Some integer identifying the success of the Evaluation
%
% Minimal example from PEtab community:
% arEvaluatePEtabSelect('modelSelectionSpecification_example_modelSelection.tsv')
%
% See also arWriteReportPEtabSelect
%
% Further information and documentation: 
% https://github.com/PEtab-dev/draft_model_selection_extension

function status = arEvaluatePEtabSelect(modelSpecFile, constrFile, config)

global ar;
arInit

%% read files
% need function like arReadPETabSelect? 

if exist('modelSpecFile','var') && ~isempty(modelSpecFile)
    modelSpecs = struct2table(tdfread(modelSpecFile));
else
    error('There is no file given that determines which models to compare.')
end

if exist('constrFile','var') && ~isempty(constrFile)
    constraintsApplied = tdfread(constrFile);
else 
    constraintsApplied = struct();
end

%% apply combinatorics of possible parameter options
% TODO, see documentation of model specification file
% basically waiting for working example that showcases this feature.

%% loop over files to be checked
status = zeros(1,height(modelSpecs));

for iLine = 1:height(modelSpecs)
    arInit

    arImportPEtab(modelSpecs.petab_yaml(iLine,:))
    ar.petab.modelID = modelSpecs.model_id(iLine,:);
    
    % apply parameter setting from row
    [pNames,idxa,idxb] = intersect(ar.pLabel,modelSpecs.Properties.VariableNames);
    if length(pNames)< width(modelSpecs)-2
        warning('Might not have gotten the full information in model specification table!')
    end
    
    for j = 1:length(pNames)
        % arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])
        if strcmp(modelSpecs.(pNames{j})(iLine,:),'estimate')
            arSetPars(pNames{j},[],1); % which ar.p at this point?
        else 
            arSetPars(pNames{j},str2double(modelSpecs.(pNames{j})(iLine,:)),0,0)
        end
    end
    
    % do fitting - LHS necessary? -- need additional fitting config?
    arFit
    
    % check for constraints
    if ~isempty(fieldnames(constraintsApplied))
        % apply conditions here.
    end
    
    % report results
    % status(iLine) = arWriteReportPEtabSelect(); % write by default all possible criteria to .yaml file
end

%% return code to check for successful run
if all (status == 0)
    status = 0;
else
    status = 1;
end

end
