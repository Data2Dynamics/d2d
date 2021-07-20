%--------------------------------------------------------------------------
% Function that calculates the objective function value in the optimization
% (as the ratio between number of model outputs and parameters, plus a 
% penalty on the number of states)
%--------------------------------------------------------------------------
function obj = objective_fun(decvector)  

global maxstates x f h p

indices   = find(decvector);
numstates = sum(decvector);

%==========================================================================
% States and functions of the reduced model:
xred = sym(zeros(numstates,1));
fred = sym(zeros(numstates,1));

for i=1:numstates
    xred(i) = x(indices(i));
    fred(i) = f(indices(i));
end

%==========================================================================
% Select outputs that measure (as many as possible of) the selected states:
hred = [];
for i=1:numel(h)
    % find the variables (states & parameters) appearing in each output:
    varsh = symvar(h(i)).';
    vars = cell(numel(varsh),1);
    for j=1:numel(varsh)
        vars{j} = char(varsh(j));
    end
    % If any of the selected states is among those variables, include that output:
    continuar = 1;
    j = 1;
    while (continuar == 1) && (j <= numel(xred))
        isstatej = strcmpi(char(xred(j)),vars);
        j = j + 1;
        if  find(isstatej) ~= 0
            hred = [hred; h(i)];
            continuar = 0;
        end
    end
end
    
%==========================================================================
% 'Variables' (parameters and states) appearing in the equations of the 
% selected states and outputs:
total_expr = [fred; hred];
variaveis = symvar(total_expr).';                       
vars = cell(numel(variaveis),1);
for i=1:numel(variaveis)
    vars{i} = char(variaveis(i));
end

% Parameters among those variables:
parslist = [];
for i=1:numel(p)
    ispari = strcmpi(char(p(i)),vars);
    if  find(ispari) ~= 0
        parslist = [parslist p(i)];
    end
end

% States (other than those in xred) among those variables:
stateslist = [];            
otherx = setdiff(x,xred); % elements of x that are not in xred        
for i=1:numel(otherx)
    ispari = strcmpi(char(otherx(i)),vars);
    if  find(ispari) ~= 0
        stateslist = [stateslist otherx(i)];
    end
end        

% The 'parameters' in the submodel include the parameters, and the
% states considered as parameters:
pred = transpose([parslist,stateslist]);
npars = numel(pred);

%==========================================================================
if numstates > maxstates
    penalty = (numstates-maxstates)^2; 
else
    penalty = 0;
end

ratio = -numel(hred)/npars;
obj = ratio + penalty;

return