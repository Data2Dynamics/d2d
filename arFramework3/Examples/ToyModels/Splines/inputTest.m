arInit;
arLoadModel('input');
arLoadData('test', 1, 'csv');
arCompileAll(true);
%arFit;
%arPlot;
title('Monotone spline (monotone between spline knots)');