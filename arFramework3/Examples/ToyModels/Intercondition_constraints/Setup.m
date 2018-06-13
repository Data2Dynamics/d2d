% This example file demonstrates how to regularize predictions in such a
% way as to minimize differences between different conditions

arInit;
arLoadModel( 'min_model' );
arLoadData( 'drug', 1, 'csv' );
arLoadData( 'noDrug', 1, 'csv'  );
arLoadData( arGenerateConditionData( 'steadystate', 'input_drug', 0, 'input_stimulus', 0 ), 1, 'csv' );
arCompileAll;

C = arFindCondition('', 'input_drug', '0', 'input_stimulus', '0' );
arSteadyState(1, C, 'all');

% Initial parameters
ar.p = -1 * ones(size(ar.p));
ar.p(arFindPar('kq')) = 3;      % Deliberately set this one wrong

