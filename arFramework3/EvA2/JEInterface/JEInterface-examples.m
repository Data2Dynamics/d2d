%% Usage examples for the EvA2 to Matlab interface.
%  Author: Marcel Kronfeld, Chair for Cognitive Systems, University of Tuebingen, Germany
%  URL: http://www.ra.cs.uni-tuebingen.de/software/EvA2/

% adapt the path settings!
addpath '/home/user/workspace/MatlabInterface' 	% .. directory containing @JEInterface
javaaddpath '/home/user/workspace/EvA2Base.jar'  	% .. the EvA2 base package
% addpath 'C:\Dokumente und Einstellungen\user\workspace\MatlabInterface'	% Windows will look differently
% javaaddpath 'C:\Dokumente und Einstellungen\user\workspace\EvA2Base.jar'	% Windows will look differently

% real valued case
R=[-5 -5 -5; 5 5 5];
JI=JEInterface(@testfun, 'double', R, R, 1, 'Display', 'iter', 'TolX', 0, 'TolFun', 0);
JI=optimize(JI, 4);
[sol, solFit]=getResult(JI);
finalPop=getMultipleResults(JI);

% binary case
R=20;
JI=JEInterface(@testfun, 'binary', R, R, 4, 'Display', 'iter');
JI=optimize(JI, 3);
[sol, fit]=getResult(JI);
finalPop=getMultipleResults(JI);

% integer case with specific initialization range
initR=[-15 -15 -15 -15 -15; -5 -5 -5 -5 -5];
R=[-15 -15 -15 -15 -15; 15 15 15 15 15];
JI=JEInterface(@testfun, 'int', R, initR, 5, 'Display', 'iter');
JI=optimize(JI, 3);
[sol, fit]=getResult(JI);
finalPop=getMultipleResults(JI);
