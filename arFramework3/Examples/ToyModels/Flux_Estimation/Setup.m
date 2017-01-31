% Example of fitting a flux to data

arInit;
arLoadModel('test');
arLoadData('test', [], 'csv');
arCompileAll(true);

arPlot;
title('Before fitting (press any key)');
ar.model.v
pause;
arFit;
arPlot;
title('After fitting');
