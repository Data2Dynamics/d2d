% Test file for splines longer than 10 knots

% The spline was generated with:
t_in=[0,1,2,4,6,9,13,20,26,70,90,110,120,130,150,160,170,180];

% Using following spline:
disp( 'Using the following spline (note that this has to be put in the def file manually):' );
splineGen(t_in)

arInit;
arLoadModel('monotone_longer');
arLoadData('test2', 1, 'csv');
ar.model(1).data(1).tExtra=ar.model(1).data(1).tExp;
arCompileAll;
arFit;
arPlot;
title('Monotone spline (monotone between spline knots)');

