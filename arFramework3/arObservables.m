%
% Observables are specified in def files in order to link the model to data
%
% Example: 
%   OBSERVABLES
%   tSTAT5_au           C   "au"  "conc."   1   1   "scale_tSTAT5 * (STAT5 + pSTAT5)"
%   pSTAT5_au           C   "au"  "conc."   1   1   "offset_pSTAT5 + scale_pSTAT5 * pSTAT5"
% 
% The arguments are specified as follows:
%   - Unique identifier (here tSTAT5_au and pSTAT5_au). Please note that 
%     this identifier should be the same as the corresponding column header 
%     in the data sheet. Note: Observable identifiers are not allowed to 
%     be the same as dynamic variable identifiers.
%   - Unit type (e.g. C for concentration).
%   - Unit itself (e.g. "au" for arbitrary units)
%   - Plain text (e.g. "conc.") for plotting purposes.
%   - Automatic scaling (0=no, 1=yes), indicates whether the maximal 
%     value in the data sheet should be rescaled to one.
%   - Logarithmic error model (0=no, 1=yes), specifies whether model and
%     data should be log10-transformed.
%   - Mathematical expression for the observable, specified in either linear 
%     or log-10 space depending on the previous argument.
% 
% In the likelihood function, the measurement noise is modeled as either a normal 
% or log-normal distribution, see 2.2 argument 6. The magnitude of the measurement 
% noise, i.e. the standard deviation of the normal distribution (or standard 
% deviation of the normal distribution in the log-space for log-normally 
% distributed noise), can be implemented as parameterized function in the
% ERRORS section.
% 
% Example:
%   ERRORS
%   tSTAT5_au		"sd_STAT5_au"
%   pSTAT5_au		"sd_STAT5_au"
%
% The arguments are specified as follows:
%   - Unique identifiers as given in the observables section (here tSTAT5_au and pSTAT5_au). 
%     It is mandatory to specify the measurement noise for each observable in the data definition file.
%
%   - Mathematical expression of the noise function (here "sd_STAT5_au" and "sd_STAT5_au" 
%     indicates that both measurements have the same measurement noise). Relative measurement 
%     noise can be implemented by using the observable identifier, e.g. 
%     "sqrt(sd_STAT5_abs^2 + sd_STAT5_rel^2 * tSTAT5_au^2)".
%
% Fixed data-based errors:
%   Sometimes, instead of using a parameterized function for the measurement 
%   noise, one would like to directly specify the amount of noise calculated 
%   from replicates in the data sheet. Put these value in the data sheet with 
%   the column header set to the observable's name but with the addition _std 
%   (here e.g. tSTAT5_au_std and pSTAT5_au_std). 
%   Please note that these values will be interpreted as standard deviation 
%   of the data not as variance! An example for this usage can be found in 
%   the example application /Example/Swameye_PNAS2003.
%
