function q_pathfound = arNEBCheckPath
% Checks the path profile after NEB calculation for connectivity

global ar

    tmpchi = ar.merger.neb.chi2_step;
    tmpresult =  ar.merger.neb.ps_result;
    
    for i = 1:size(tmpresult,1)-1
        tmpstepsize(i) = norm(tmpresult(i,:) - tmpresult(i+1,:));
    end

    threshold_tolerance = max( abs(tmpchi(1)-tmpchi(end))/100, 1e-4);
    
    % below start and end?   
    q_below_chi2 = sum(tmpchi > ...
        max(tmpchi(1),tmpchi(end))) < 1;
    
    % tolerance for monoticitiy test
    q_monotonic = check_monotonicity_path(tmpchi, ...
        threshold_tolerance);
    
    % stepsize not too extreme
    q_stepextreme = ...
    sum( tmpstepsize > 0.3 * ( norm(tmpresult(1,:) - tmpresult(end,:))) )  < 1;
    
    q_pathfound = sum(q_stepextreme * q_monotonic * q_below_chi2) > 0;
    

end

