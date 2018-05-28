% Example of fitting a flux to data
Setup

%%
arPlot;
title('Before fitting (press any key)');
ar.model.v
pause;
arFit;
arPlot;
title('After fitting');
