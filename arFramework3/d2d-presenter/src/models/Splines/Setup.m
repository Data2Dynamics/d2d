% Test file for different spline functions

arInit;
arLoadModel('normal_cubic');
arLoadData('test', 1, 'csv');
arLoadModel('positive_cubic');
arLoadData('test', 2, 'csv');
arLoadModel('monotone');
arLoadData('test', 3, 'csv');
arCompileAll;
arFit;
%arPlot;
%title('Cubic spline');
%pause;