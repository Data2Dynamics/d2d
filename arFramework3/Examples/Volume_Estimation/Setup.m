% Example of fitting cell volumes
%
% Note that currently, this resorts to using FD for the inner sensitivities
% (sensitivies used within the ODE solver), but sensitivity equations for
% the outer ones.

arInit;
arLoadModel('test');
arLoadData('test', [], 'csv');
arCompileAll(true);

arPlot;
title('Before fitting');

ar.qFit=ones(size(ar.qFit));
arFit;

arPlot;
title('After fitting');
