function int = postProcess(int, steps, sigmaClust, varargin)
% Do post processing of the last optimization results and export
% the solution set to the JEInterface.
%   int=postProcess(int, steps, sigmaClust [, nBest])
% Either of the parameters steps and sigmaClust may be 0 to leave out
% hill-climbing or clustering step. Setting both to 0 is not meaningful.
% The optional nBest parameter gives an upper bound to the number of
% returned solutions, which are sorted by fitness. Leaving this parameter
% out means that all solutions suggested by the optimizer (after
% clustering if activated) are returned, the number of which is usually,
% but not for all optimizers (e.g., not for CBN), limited to the population size.
%
% Arguments:
%   int: the JEInterface instance
%   steps: number of hill climbing steps to perform or 0
%   sigmaClust: paramter for the density based clustering to perform or 0,
%       relative to the problem range
%   nBest: (optional) maximum number of solutions to return


%   stepSize: (optional) step size of the stochastic hill climber. 

if (int.finished == 0)
    error('please wait for the current run to finish');
end

int.finished = 0;

if (nargin > 3) && (isnumeric(varargin{1}))
    nBest = varargin{1};
else
	nBest = -1;
end

int=runEvalLoopJE(int, 2, -1, '', steps, sigmaClust, nBest);
