%PEtab Problem class

classdef Problem < handle & matlab.mixin.Copyable & dynamicprops
    %PEtab parameter estimation problem as defined by:
    %
    %- SBML model
    %- condition table
    %- measurement table
    %- parameter table
    %- observables table
    %
    %Optionally it may contain visualization tables.
    %
    %Properties:
    %   sbml_model [Sbml object]:
    %       SBML model.
    %   condition_df [table]: 
    %       PEtab condition table.
    %   measurement_df [table]: 
    %       PEtab measurement table.
    %   parameter_df [table]: 
    %       PEtab parameter table.
    %   observable_df [table]: 
    %       PEtab observables table.
    %
    %   amimodel [struct]:
    %       AMICI model struct.
    
    properties (SetAccess = private)
        folder string = [];
    end
    
    properties
        sbml_model Sbml = Sbml.empty();
        condition_df table = table.empty();
        measurement_df table = table.empty();
        parameter_df table = table.empty();
        observable_df table = table.empty();
        
        amimodel struct = struct.empty();
    end
    
    methods
        %% CLASS CONSTRUCTOR
        function obj = Problem(folder, sbml_model, condition_df, ...
                measurement_df, parameter_df, observable_df) %-> Problem
            %Problem constructor.
            
            if nargin == 0
                return
            end
            
            obj.folder = folder;
            
            obj.sbml_model = sbml_model;
            obj.condition_df = condition_df;
            obj.measurement_df = measurement_df;
            obj.parameter_df = parameter_df;
            obj.observable_df = observable_df;
        end
        
        %% OBJECT METHODS        
        function getDynamics(obj)
            %Loads AMICI model structure containing model dynamics.
            
            fprintf('Loading dynamics...\n')
            
            tmp = SBMLode(obj.sbml_model.sbml);
            tmp.writeAMICI(obj.sbml_model.id, obj.folder);            
                        
            obj.amimodel = feval([obj.sbml_model.id ...
                '_syms']);
            
            filename = fullfile(obj.folder, [obj.sbml_model.id ...
                '_syms.m']);            
            delete(filename);
            
            fprintf('Dynamics succesfully loaded.\n\n')
        end
        
        function out = get_optimization_parameters(obj)
            %Return string array of optimization parameter IDs.
            %
            %See parameter/get_optimization_parameters
            
            out = get_optimization_parameters(obj.parameter_df);
        end
        
        function out = get_optimization_parameter_scales(obj)
            %Return dictionary of optimization parameter scales.
            %
            %See parameters/get_optimization_parameter_scales
            
            out = get_optimization_parameter_scales(obj.parameter_df);
        end
        
        function out = get_model_parameters(obj)
            %See sbml/get_model_parameters
            
            out = get_model_parameters(obj.sbml_model);
        end
        
        function out = get_observable_ids(obj)
            %Returns string array of observable ids.
            
            out = transpose(string(obj.observable_df.observableId));
        end
        
        function out = x_ids(obj) %-> [string]
            %Parameter table parameter IDs.
            
            out = obj.get_x_ids;
        end
        
        function out = x_free_ids(obj) %-> [string]
            %Parameter table parameter IDs, for free parameters.
            
            out = obj.get_x_ids(true, false);
        end
        
        function out = x_fixed_ids(obj) %-> [string]
            %Parameter table parameter IDs, for fixed parameters.
            
            out = obj.get_x_ids(false);
        end
        
        function out = x_nominal(obj)
            %Parameter table nominal values.
            
            out = obj.get_x_nominal;
        end
        
        function out = x_nominal_free(obj)
            %Parameter table nominal values, for free parameters.
            
            out = obj.get_x_nominal(true, false);
        end
        
        function out = x_nominal_fixed(obj)
            %Parameter table nominal values, for fixed parameters.
            
            out = obj.get_x_nominal(false, true);
        end
        
        function out = x_nominal_scaled(obj)
            %Parameter table nominal values with applied parameter scaling.
            
            out = obj.get_x_nominal(true, true, true);
        end
        
        function out = x_nominal_free_scaled(obj)
            %Parameter table nominal values with applied parameter scaling,
            %for free parameters.
            
            out = obj.get_x_nominal(true, false, true);
        end
        
        function out = x_nominal_fixed_scaled(obj)
            %Parameter table nominal values with applied parameter scaling,
            %for fixed parameters.
            
            out = obj.get_x_nominal(false, true, true);
        end
        
        function out = lb(obj)
            %Parameter table lower bounds.
            
            out = obj.get_lb;
        end
        
        function out = lb_scaled(obj)
            %Parameter table lower bounds with applied parameter scaling.
            
            out = obj.get_lb(true, true, true);
        end
        
                function out = ub(obj)
            %Parameter table lower bounds.
            
            out = obj.get_ub;
        end
        
        function out = ub_scaled(obj)
            %Parameter table lower bounds with applied parameter scaling.
            
            out = obj.get_ub(true, true, true);
        end
        
        function out = get_simulation_conditions_from_measurement_df(obj)
            %See measurements/get_simulation_conditions
            
            out = get_simulation_conditions(obj.measurement_df);
        end
        
        function [preeq, sim] = ...
                get_optimization_to_simulation_parameter_mapping(obj, ...
                warn_unmapped, scaled_parameters)
            %See parameter_mapping/get_simulation_to_optimization_parameter_mapping
            
            if nargin == 2
                scaled_parameters = false;
            elseif nargin == 1
                warn_unmapped = true;
                scaled_parameters = false;
            end
            
            [preeq, sim] = ...
                get_optimization_to_simulation_parameter_mapping( ...
                obj.condition_df, obj.measurement_df, 'parameter_df', ...
                obj.parameter_df, 'observable_df', obj.observable_df, ...
                'sbml_model', obj.sbml_model, 'warn_unmapped', ...
                warn_unmapped, 'scaled_parameters', scaled_parameters);
        end
    end
    
    methods (Static)
        %% STATIC METHODS
        function out = from_files(varargin) %-> [Problem]
            %Loads model and tables from files.
            %
            %Arguments:
            %   folder Optional string/*["."]:
            %       Model folder.
            %   sbml_file Optional string/*[]:
            %       SBML model.
            %   condition_file Optional string/*[]: 
            %       PEtab condition table.
            %   measurement_file Optional string/*[]: 
            %       PEtab measurement table.
            %   parameter_file Optional string/*[]: 
            %       PEtab parameter table.
            %   observables_file Optional string/*[]: 
            %       PEtab observables table.
            %
            %Returns:
            %   [Problem]
            %       PEtab Problem object.
            
            p = inputParser;
            addParameter(p, 'folder', '.')
            addParameter(p, 'sbml_file', Sbml.empty())
            addParameter(p, 'condition_file', table.empty())
            addParameter(p, 'measurement_file', table.empty())
            addParameter(p, 'parameter_file', table.empty())
            addParameter(p, 'observable_file', table.empty())
            parse(p, varargin{:})
            
            folder = what(p.Results.folder).path;
            
            sbml_file = p.Results.sbml_file;
            condition_file = p.Results.condition_file;
            measurement_file = p.Results.measurement_file;
            parameter_file = p.Results.parameter_file;
            observable_file = p.Results.observable_file;            
            
            if ~isempty_ext(sbml_file)
                sbml_model = Sbml(sbml_file);
            end
            
            if ~isempty_ext(condition_file)
                condition_df = get_condition_df(condition_file);
            end
            
            if ~isempty(measurement_file)
                measurement_df = concat_tables(measurement_file, ...
                    @get_measurement_df);
            end
            
            if ~isempty(parameter_file)
                parameter_df = get_parameter_df(parameter_file);
            end
            
            if ~isempty(observable_file)
                observable_df = concat_tables(observable_file, ...
                    @get_observable_df);
            end
            
            out = Problem(folder, sbml_model, condition_df, ...
                measurement_df, parameter_df, observable_df);
        end
        
        function out = from_yaml(folderOrPath, modelname) %-> [Problem]
            %Load model and tables as specified by YAML file.
            %
            %Arguments:
            %   folder string:
            %       PEtab configuration as dictionary or Path to the 
            %       directory in which the yaml file is located.
            %   modelname Optional string/*[]:
            %       If specified, overrides the model component in the file
            %       names. Defaults to the last component of 'folder'.
            %
            %Returns:
            %   [Problem]
            %       PEtab Problem object.            
            
            if nargin == 1
                modelname = [];
            end
            
            folderOrPath = char(folderOrPath);
            modelname = char(modelname);
            
            if ~exist(folderOrPath)
                error('FROM_FOLDER:FileNotFoundError', ...
                    'No such file or directory')
            end
            
            if isempty_ext(modelname)
                [folder, ~] = fileparts(folderOrPath);
            else
                folder = folderOrPath;
                folderOrPath = fullfile(folder, modelname + ".yaml");
            end
            
            yaml_config = ReadYaml(char(folderOrPath));
            
            problem0 = yaml_config.problems{1};
            
            assert_single_condition_and_sbml_file(problem0);
            
            out = Problem.from_files('folder', folder, ...
                'sbml_file', ...
                fullfile(folder, problem0.sbml_files{1}), ...
                'condition_file', ...
                fullfile(folder, problem0.condition_files{1}), ...
                'measurement_file', ...
                fullfile(folder, problem0.measurement_files{1}), ...
                'parameter_file', ...
                fullfile(folder, yaml_config.parameter_file), ...
                'observable_file', ...
                fullfile(folder, problem0.observable_files{1}));
        end
    end
        
    methods (Access = public)
        %% AUXILIAR METHODS
        function out = x_free_indices(obj)
            %Parameter table estimated parameter indices.
            
            out = transpose(find(obj.parameter_df.estimate));
        end
        
        function out = x_fixed_indices(obj)
            %Parameter table non-estimated parameter indices.
            
            out = transpose(find(~obj.parameter_df.estimate));
        end       
        
        function out = apply_mask(obj, v, free, fixed) %-> [[string], [cell]]
            %Apply mask of only free or only fixed values.
            %
            %Arguments:
            %   v [[string], [cell]]:
            %       The full array the mask is to be applied to.
            %   free bool:
            %       Whether to return free parameters, i.e., parameters to
            %       estimate.
            %   fixed bool:
            %       Whether to return fixed parameters, i.e. parameters not
            %       to estimate.
            %
            %Returns:
            %   [[string], [cell]]
            %       The reduced vector with applied mask.
            
            if nargin == 2
                fixed = true;
            elseif nargin == 1                
                free = true;
                fixed = true;
            end
            
            if ~free && ~fixed
                out = [];
            elseif ~free
                out = v(obj.x_fixed_indices);
            elseif ~fixed
                out = v(obj.x_free_indices);
            else
                out = v;
            end
        end
        
        function out = get_x_ids(obj, free, fixed) %-> [string]
            %Generic function to get parameter ids.
            %
            %Parameters:
            %   free bool:
            %       Whether to return free parameters, i.e. parameters to 
            %       estimate.
            %   fixed bool:
            %       Whether to return fixed parameters, i.e. parameters 
            %       not to estimate.
            %
            %Returns
            %   [string]
            %       The parameter ids.
            
            if nargin == 2
                fixed = true;
            elseif nargin == 1
                free = true;
                fixed = true;
            end
            
            v = transpose(string(obj.parameter_df.parameterId));
            out = obj.apply_mask(v, free, fixed);
        end
        
        function out = get_x_nominal(obj, free, fixed, scaled) %-> [numeric]
            %Generic function to get parameter nominal values.
            %
            %Parameters:
            %   free bool:
            %       Whether to return free parameters, i.e. parameters to 
            %       estimate.
            %   fixed bool:
            %       Whether to return fixed parameters, i.e. parameters 
            %       not to estimate.
            %   scaled bool:
            %       Wheter to scale the values according to the parameter
            %       scale, or return them on linear scale.
            %
            %Returns
            %   [numeric]
            %       The parameter nominal values.
            
            if nargin == 3
                scaled = false;
            elseif nargin == 2
                fixed = true;
                scaled = false;
            elseif nargin == 1
                free = true;
                fixed = true;
                scaled = false;
            end
            
            v = obj.parameter_df.nominalValue;
            
            if scaled                
                v = map(@scale, num2cell(v), ...
                    obj.parameter_df.parameterScale);
            end
            
            out = transpose(obj.apply_mask(v, free, fixed));
        end
        
        function out = get_lb(obj, free, fixed, scaled) %-> [numeric]
            %Generic function to get lower parameter bounds.
            %
            %Parameters:
            %   free bool:
            %       Whether to return free parameters, i.e. parameters to 
            %       estimate.
            %   fixed bool:
            %       Whether to return fixed parameters, i.e. parameters 
            %       not to estimate.
            %   scaled bool:
            %       Wheter to scale the values according to the parameter
            %       scale, or return them on linear scale.
            %
            %Returns
            %   [numeric]
            %       The lower parameter bounds.
            
            if nargin == 3
                scaled = false;
            elseif nargin == 2
                fixed = true;
                scaled = false;
            elseif nargin == 1
                free = true;
                fixed = true;
                scaled = false;
            end
            
            v = obj.parameter_df.lowerBound;
            
            if scaled                
                v = map(@scale, num2cell(v), ...
                    obj.parameter_df.parameterScale);
            end
            
            out = transpose(obj.apply_mask(v, free, fixed));
        end
        
        function out = get_ub(obj, free, fixed, scaled) %-> [numeric]
            %Generic function to get upper parameter bounds.
            %
            %Parameters:
            %   free bool:
            %       Whether to return free parameters, i.e. parameters to
            %       estimate.
            %   fixed bool:
            %       Whether to return fixed parameters, i.e. parameters
            %       not to estimate.
            %   scaled bool:
            %       Wheter to scale the values according to the parameter
            %       scale, or return them on linear scale.
            %
            %Returns
            %   [numeric]
            %       The upper parameter bounds.
            
            if nargin == 3
                scaled = false;
            elseif nargin == 2
                fixed = true;
                scaled = false;
            elseif nargin == 1
                free = true;
                fixed = true;
                scaled = false;
            end
            
            v = obj.parameter_df.upperBound;
            
            if scaled
                v = map(@scale, num2cell(v), ...
                    obj.parameter_df.parameterScale);
            end
            
            out = transpose(obj.apply_mask(v, free, fixed));
        end
    end
end