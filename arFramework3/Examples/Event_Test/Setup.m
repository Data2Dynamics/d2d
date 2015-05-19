% Load models & data
arInit;
arLoadModel('pulseTest');
arLoadData('pulse',1,'csv',true);

arCompileAll;

% Plot the results without events
arPlot; title( 'Without events [press any key]' );
disp('Press any key' );

pause;

% Find events belonging to step functions and set them as reinitialization
% time points.
arFindInputs;

% Plot the results with events
arPlot; title( 'With events' );
