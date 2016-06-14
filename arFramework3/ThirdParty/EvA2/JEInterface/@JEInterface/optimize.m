function retInt = optimize(int, optType, varargin)
% Start a EvA2 optimization run.
%       optimize(interface, optType, [, outputFilePrefix ] )
% where
%       interface: instance of JEInterface
%       optType: integer indicating the type of the optimization strategy
%       to use.
%       resultFilePrefix: (optional) char prefix for an optional verbose
%           output file

if (int.finished == 0) 
    error('please wait for the current run to finish');
end
if ((nargin == 2) || (nargin == 3))
    if (nargin == 3) 
        outputFilePrefix = varargin{1};
    else
        outputFilePrefix = 'none';
    end

    if (~isa(int, 'JEInterface') || ~isscalar(optType) || ~isa(outputFilePrefix, 'char'))
        error('Invalid argument!')
    end
    int.finished = 0;
    int.msg = 'running...';
    % adapt options possibly changed by the user, concerning
    xTol = int.opts.TolX;
    maxEvals = int.opts.MaxFunEvals;
    fTol = int.opts.TolFun;
    if (ischar(xTol)) ;     xTol    =   str2num(xTol); end;
    if (ischar(maxEvals)) ; maxEvals=   str2num(maxEvals); end;
    if (ischar(fTol)) ;     fTol    =   str2num(fTol); end;

    import eva2.server.go.operators.terminators.PhenotypeConvergenceTerminator;
    import eva2.server.go.operators.terminators.FitnessConvergenceTerminator;
    import eva2.server.go.operators.terminators.PopulationMeasureTerminator;  
    import eva2.server.go.operators.terminators.PopulationMeasureTerminator.*;
    import eva2.server.go.operators.terminators.CombinedTerminator;
    import eva2.server.go.operators.terminators.EvaluationTerminator;
    import eva2.OptimizerFactory;
    import eva2.server.go.problems.MatlabProblem;

    % set some default values if theyre not given
    % fminsearch, for example, always uses TolX and TolFun with default
    % values of 1e-4 in . Thats what we do as well
    if (isempty(int.opts.TolX)) ; xTol = 1e-4; end
    if (isempty(int.opts.TolFun)) ; fTol = 1e-4; end
    
    % construct Terminators
    if ((xTol > 0) && (fTol > 0))
        % both criteria are given, use combination
        convTerm = CombinedTerminator(MatlabProblem.makeFitConvTerm(fTol, int.opts.TolFunEvals), ...
                    MatlabProblem.makePhenConvTerm(xTol, int.opts.TolXEvals), 1);
    else if (xTol > 0)  % only phenotye convergence
            convTerm = MatlabProblem.makePhenConvTerm(xTol, int.opts.TolXEvals);
        else if (fTol > 0 )             % only fitness covnergence
                convTerm = MatlabProblem.makeFitConvTerm(fTol, int.opts.TolFunEvals);
            else
               convTerm = 'undef'; % signal that there is no terminator yet
            end
        end
    end

    if (ischar(convTerm)) % if no convergence terminator is defined so far, use fitness calls 
        if (isempty(maxEvals))
            error('Error: no termination criterion defined! Please check options.');
            % int.opts.MaxFunEvals = OptimizerFactory.getDefaultFitCalls;
            % maxEvals = OptimizerFactory.getDefaultFitCalls;
        end
        convTerm = EvaluationTerminator(maxEvals);
        eva2.OptimizerFactory.setTerminator(convTerm);
    else % there is a convergence terminator              
        eva2.OptimizerFactory.setTerminator(convTerm); % so set it
        if (~isempty(maxEvals) && (maxEvals > 0))
            % if TolX/TolFun plus MaxFunEvals is defined additionally, combine an
            % EvaluationTerminator in disjunction, as Matlab does.
            eva2.OptimizerFactory.addTerminator(EvaluationTerminator(maxEvals), 0);
        end
    end
    int.mp.setOutputAllStatFields(int.outputAllStatsFields==1);

    % set display
    if (strcmp(int.opts.Display,'off') || isempty(int.opts.Display))
        int.mp.setStatsOutput(0);
    elseif (strcmp(int.opts.Display, 'final'))
        int.mp.setStatsOutput(1);
    elseif (strcmp(int.opts.Display, 'notify')) 
        % 'notify' is not the perfect notion for "show every k-th
        % iteration", but optimset wont allow changing names.
        int.mp.setStatsOutput(2);       
    elseif (strcmp(int.opts.Display, 'iter'))
        % this should rather be 2 in JE slang, but in matlab slang its more like 3
        int.mp.setStatsOutput(3);    
    else
        error('invalid Display option, only off/final/notify/iter are recognized');
    end
    
    if isempty(int.seedPop) % set the seed data
        int.mp.clearSeedPopulation;
    else
        int.mp.setSeedPopulation(int.seedPop, int.seedPopFit);
    end
    int=runEvalLoopJE(int, 1, optType, outputFilePrefix, -1, -1, -1);
    
else
    error('Wrong number of arguments!')
end
retInt=int;
