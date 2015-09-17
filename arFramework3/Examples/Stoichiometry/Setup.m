% Small test file for reactions with non unity stoichiometry

% Various test examples with CUSTOM
arInit;
arLoadModel('test');
arCompileAll;
arPlot;

pause;

% Comparing MASSACTION with CUSTOM
arInit;
arLoadModel('test2');
arCompileAll;
arPlot;