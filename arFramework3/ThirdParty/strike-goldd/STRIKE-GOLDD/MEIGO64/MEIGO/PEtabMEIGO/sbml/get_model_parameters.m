
function out = get_model_parameters(sbml_model, with_values) %-> [string, [Dict]]
    %Return SBML model parameters which are not AssignmentRule targets for
    %observables or sigmas.
    %
    %Arguments:
    %   sbml_model [libsbml struct]: 
    %       SBML model struct.
    %   with_values Optional bool/*false: 
    %       If false, returns list of SBML model parameter IDs which are 
    %       not AssignmentRule targets for observables or sigmas. If true,
    %       returns a dictionary with those parameter IDs as key and
    %       parameter values from the SBML model as values.
    %
    %Returns:
    %   [string/[Dict]]
    %       List of model parameters ids or dictionary with model 
    %       parameter ids/values pairs.
    
    if nargin == 1
        with_values = false;
    end
    
    testfunc = @(id) isempty_ext( ...
        sbml_model.getAssignmentRuleByVariable(id));    
    mask = map(testfunc, {sbml_model.parameter.id});
    
    ids = string({sbml_model.parameter(mask).id});
    vals = {sbml_model.parameter(mask).value};
    
    if with_values
        out = Dict(ids, vals);
    else
        out = ids;
    end
end