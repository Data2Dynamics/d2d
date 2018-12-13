%% Check
ar.config.atol=1e-8;
ar.config.rtol=1e-8;
arQplot('x')
arPlot

arCompareWithBiobaseSimulation('SIMU1447954318090.dat');
