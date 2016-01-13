% Test file for different spline functions

arInit;
arLoadModel('normal_cubic');
arLoadData('test', 1, 'csv');
arCompileAll;
arFit;
arPlot;
title('Cubic spline');
pause;

arInit;
arLoadModel('positive_cubic');
arLoadData('test', 1, 'csv');
arCompileAll;
arFit;
arPlot;
title('Cubic spline (forced positive)');
pause;

arInit;
arLoadModel('monotone');
arLoadData('test', 1, 'csv');
arCompileAll;
arFit;
arPlot;
title('Monotone spline (monotone between spline knots)');