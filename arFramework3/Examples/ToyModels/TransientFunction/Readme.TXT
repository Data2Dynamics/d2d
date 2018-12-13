The "transient function" is a universal tool describing many time course from the systems biology field data reasonably. This folder illustrates how the "transient function" can be fitted to exp. data.

Fitting involves testing both directions: activatory and inhibitory dynamics. For fitting both possibilitites is done by arFitTransient. For using arFitLHS.m, arFits had to be adapted. Moreover, reasonable bounds are required especially for the time-scales. This is done in DefaultLbUbTransient.m . The meaning of the parameters is determined in Initialize_FitTransient.m. All these functions are in folder TransientFunction_library.

The examples comprise the following settings:
1) Example1_FitFunctionToData: 
Fitting several experiments with individual parameters. Also data of untreated conditions (assuming a constant dynamics) might be included.

