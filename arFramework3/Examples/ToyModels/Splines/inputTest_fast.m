% This file shows an example of how to use fixed input splines. These
% splines cannot take free parameters, but they can be used to import
% input functions with a large number of time points.

arInit;
ar.config.turboSplines = 1; % Enable spline coefficient caching
arLoadModel('input');
arLoadData('long', 1, 'csv');
arCompileAll(true);

% Make sure that we actually see every data point (otherwise it seems that
% not all the data is fitted, because the displayed points do not end up at
% the spline knots.
ar.model.data.tExtra = [1:5000]; arLink;

% We should really use the new plots because the performance of the old
% ones is terrible when having so many points.
ar.config.useNewPlots=1;
arPlot;
title('Monotone spline (monotone between spline knots)');