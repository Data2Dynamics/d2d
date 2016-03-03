% Small test file for reactions with non unity stoichiometry

findState = @(ar, state)find( strcmp(ar.model.x, state ) );

% Various test examples with CUSTOM
arInit;
arLoadModel('test');
arCompileAll;
%arPlot;

X = ar.model.condition.xFineSimu(:,end);
%pause;

% Comparing MASSACTION with CUSTOM
arInit;
arLoadModel('test2');
arCompileAll;
%arPlot;
