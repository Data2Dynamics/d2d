function [default]=vns_defaults

%Assings default values for all the options
default.maxeval                   =       1000;            %Maximum number of function evaluations
default.maxtime                   =       60;              %Maximum CPU time
default.maxdist                   =       0.5;             %Percentage of the problem dimension which will be perturbed in the furthest neighborhood (vary between 0-1)
default.use_local                 =       1;               %Uses local search (1) or not (0)

%The following options only apply if use_local is active
default.aggr                      =       0;               %Aggressive search. The local search is only applied when the best solution has been improved (1=aggressive search, 0=non-aggressive search)
default.local_search_type         =       1;               %Applies a first (=1) or a best (=2) improvement scheme for the local search
default.decomp                    =       1;               %Decompose the local search (=1) using only the variables perturbed in the global phase



