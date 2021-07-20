function apply_parameter_table(par_mapping, scale_mapping, ...
    parameter_df, scaled_parameters)
    % Replace parameters from parameter table in mapping list for a given
    % condition and set the corresponding scale.
    %
    % Replace non-estimated parameters by "nominalValues"
    % (un-scaled / lin-scaled), replace estimated parameters by the 
    % respective ID.
    %
    % Parameters:
    %   par_mapping (Dict):
    %        dict obtained from "get_parameter_mapping_for_condition".
    %   scale_mapping (Dict):
    %        Parameter scales dict obtained from 
    %        "get_parameter_mapping_for_condition".
    %   parameter_df (table):
    %       PEtab parameter table.
    %   scaled_parameters (bool):
    %       Whether to scale parameters or not.
    
    if nargin == 3
        scaled_parameters = false;        
    elseif nargin == 2
        parameter_df = table.empty();
        scaled_parameters = false;
    end
    
    if isempty_ext(parameter_df)
        return
    end
    
    parids = transpose(intersect(par_mapping.keys, ...
        parameter_df.parameterId, 'stable'));
    parameter_df.Properties.RowNames = parameter_df.parameterId;
    parameter_df.parameterId = [];    
    
    scales = transpose(parameter_df.parameterScale(parids));    
    scale_mapping.addpairs(parids, cellstr(scales));
    
    estpars = transpose(parameter_df.Properties.RowNames( ...
        logical(parameter_df.estimate)));
    nonestpars = transpose(setdiff(parameter_df.Properties.RowNames, ...
        estpars, 'stable'));
    
    vals = transpose(parameter_df.nominalValue(nonestpars));
    if scaled_parameters
        scales = transpose(parameter_df.parameterScale(nonestpars));
        vals = map(@scale, num2cell(vals), scales);
    else
        scale_mapping.addpairs(parids, cellstr(repmat("lin", 1, ...
            numel(parids))));
    end
    
    par_mapping.addpairs(nonestpars, vals);
    par_mapping.addpairs(string(estpars), estpars);
    
    for i = 1:numel(par_mapping.keys)
        problem_par = par_mapping.keys{i};
        sim_par = par_mapping.values{i};
        
        if isnumeric(sim_par)
            continue
        end
        
        try            
            par_mapping(problem_par) = par_mapping(sim_par);
            scale_mapping(problem_par) = scale_mapping(sim_par);
        catch
            if isempty_ext(parameter_df)
                error([]);
            end
            
            if ismember('parameterScale', ...
                    parameter_df.Properties.VariableNames)
                scalestr = parameter_df.parameterScale{sim_par};
            else
                scalestr = 'lin';
            end
            
            if ismember('estimate', ...
                    parameter_df.Properties.VariableNames) && ...
                        ~parameter_df.estimate(sim_par)
                
                val = parameter_df.nominalValue(sim_par);
                if scaled_parameters
                    val = scale(val, scalestr);
                else
                    scalestr = 'lin';
                end
                
                par_mapping(problem_par) = val;
            end
            
            scale_mapping(problem_par) = scalestr;
        end
    end
end