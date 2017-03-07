% This file tests optimization of models with complex step functions
% It is expected to throw a warning about the location parameter. Since
% we are not optimizing over the location parameter, this is ok though.
%
% See estimateLocaton.m for an example where the location is also unknown.

arInit;
arLoadModel('stepInput');
arLoadData('stepInput',1,'csv',true);

%% compile
arCompileAll();

ar.config.useEvents = 1;

arSimu(true,true,true);
arGetMerit();

