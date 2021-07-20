function out = getMeigoProblem(folder, modelname, scaled_bounds)
    %TODO
    
    out = struct();
    
    if nargin == 1
        petab = Problem.from_yaml(folder);
        scaled_bounds = false;
    elseif nargin == 2
        petab = Problem.from_yaml(folder, modelname);
        scaled_bounds = false;
    end
    
    out.f = @(p) calculateProblemDesviation(petab, p);
    
    if scaled_bounds
        out.x_L = petab.lb_scaled;
        out.x_U = petab.ub_scaled;        
    else
        out.x_L = petab.lb;
        out.x_U = petab.ub;
    end
    
    out.x_0 = petab.x_nominal_free;    
end