function [sol, fit] = getResult(int)
% Returns the optimization solution and its fitness if the run has been finished, or an
% intermediate solution and ints fitness if the run has not finished yet or an empty array
% if there is no intermediate solution yet.
%       [sol, fit] = getResult(int)

if (isFinished(int)) 
    sol = int.result';
else
    sol = int.mp.getIntermediateResult()';
end

if (isempty(sol)) 
    fit = NaN;
else
    if (isempty(int.range)) 
        sol=convertUnsignedJE(int, sol);
    end;
        
    if (isempty(int.args))
        fit = feval(int.f, sol);
    else
        fit = feval(int.f, sol, int.args);
    end
end
%disp('Fitness: ');
%disp(y);