The "transient function" is a universal tool describing many time course from the systems biology field data reasonably. This folder illustrates how the "transient function" can be fitted to exp. data.

The data in this example consists of three experiments, each experiment has two time courses (treated and untreated). In this example, individual parameters were created for each experiment. Each experiment has the following parameters: 2 amplitudes (for sustained and transient parts), 2 time-scales, 1 offset and 1 error-parameter.

If you want to fit your own data, just replace the example data. Leave the column names nExpID	isTreated unchanged. If you only have a single experiment or treatment, then the respective columns just have the same value in every row. If you have more treatments (than "treated" and "untreated"), this should be handled by defining extra observation functions. Sorry, I dont't want to explain this further at this point. 

If "untreated" data is not available, then the corresponding rows in the data table are missing. The code/example should work anyway.

Implementation of a common error for all experiments is illustrated in ExampleData.def.

Fitting involves testing both directions: activatory and inhibitory dynamics. For fitting both possibilitites is done by arFitTransient. For using arFitLHS.m, arFits had to be adapted. Moreover, reasonable bounds are required especially for the time-scales. This is done in DefaultLbUbTransient.m . The meaning of the parameters is determined in Initialize_FitTransient.m. All these functions are in folder TransientFunction_library.

In case of troubble, contact Clemens Kreutz, ckreutz@fdm.uni-freiburg.de