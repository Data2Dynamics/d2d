function int = setOutputAllStatsFields(int, doOutputAllFields)
% (De-)Activate printing all available statistic data to the text log. Deactivation may
%	improve performance for larger population sizes.
%       int = setOutputAllStatsFields(int, doOutputAllFields)
%       	int: instance of JEInterface
%       	doOutputAllFields: 1 or 0 for activation or deactivation, respectively

int.outputAllStatsFields=doOutputAllFields;
