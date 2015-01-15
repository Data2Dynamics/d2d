function int=runEvalLoopJE(int, optOrPostProc, optType, outputFilePrefix, steps, sigmaClust, nBest)
% Internal method starting a EvA2 optimization loop.
% Calling this directly may interfere with optimization.

% This function handles the communciation between JE and Matlab, main
% optimization configuration is done in the
% optimize/optimizeWith/postProcess functions from which this one should be
% called.
% optOrPostProc: 1 for optimize, 2 for postProcess
% optType, outputFilePrefix are parameters for optimize, dont care when
% postprocessing.
% steps, sigmaClust and nBest are parameters for postProcess, dont care
% when optimizing. nBest may be -1 to show all.

global stopOptimization
global JEMediator

if ~isempty(int.mediator)
    int.mediator.quit
    int.mediator='';
end
% disp(sprintf('creating mediator'));

% set up a mediator and inform JE
int.mediator = eva2.server.go.problems.MatlabEvalMediator(int.opts.NiceSleepTime);
int.mp.setMediator(int.mediator);
JEMediator=int.mediator;
createStopBox=int.opts.CreateStopBox;

% disp(sprintf('mediator created, calling optimize'));

% start the JE thread
if (optOrPostProc == 1)
    stopText='Stop EvA2 optimization';
    int.mp.optimize(optType, outputFilePrefix, int.optParams, int.optParamValues);
else % post processing
    stopText='Stop EvA2 post processing';
    int.mp.requestPostProcessing(steps, sigmaClust, nBest);
end
% disp(sprintf('calling optimize done'));

% handle the case when the optimization has been called from a users script
% and not from the toolboxes parameter estimation function (which has an
% own stop button). we decide this by checking if the global variable
% stopOptimization is empty. if it is then it is not the toolbox calling
% and we create an own button to stop it.
if isempty(stopOptimization),    
    % set switch to 0 now
    stopOptimization = 0;
    startTime=clock;
    timeStr=sprintf('%d.%d. %d:%0.2d:%0.2d',startTime(3), startTime(2), startTime(4), startTime(5), round(startTime(6)));
    disp(sprintf('Starting optimization at %s', timeStr));
    % create a cancel button box (case without SBtoolbox)
    stopText=sprintf('%s (%s)', stopText, timeStr);
    if (createStopBox == 1)
        boxHandle=figure('Position',[100 600 250 80], 'MenuBar', 'none', 'Name', 'EvA2 optimization...', 'NumberTitle','off');
        uicontrol(boxHandle,'Style', 'pushbutton', 'String', 'Cancel', 'Position', [25 25 60 30], 'Callback', 'global stopOptimization; stopOptimization=1;');
        uicontrol(boxHandle,'Style', 'text', 'String', stopText, 'Position', [100 15 130 50]);
        drawnow;
%    else
%        disp(stopText);
    end
   
    % set flag for non toolbox optimization
    nontoolboxopt = 1;
else
    % disp('seems like the toolbox is going on');
    % its an estimation using the toolbox' parameter estimation thing
    nontoolboxopt = 0;
end

stopOnce=1;
cnt=1;
% disp(sprintf('before eval loop... %d',cnt));
% repeat the mediator thread and eval call until finished
try
    while (~int.mediator.isFinished())
 		% disp(sprintf('running mediator id %d',cnt));
        int.mediator.run(cnt);
 		% disp(sprintf('after running mediator id %d',cnt));
        cnt=cnt+1;
        if (~int.mediator.isFinished())
		% disp('getting question');
            x = int.mediator.getQuestion();
		%disp('question asked');
            if (isempty(int.range))
                %size(x)
                x=convertUnsignedJE(int, x);
                %disp('here B');
                %x
            end
%            size(x)
            try
                if (isempty(int.args))
                    res = feval(int.f, x);
                else
                    res = feval(int.f, x, int.args);
                end
 		% disp(sprintf('res is %d',res));
                %res
            catch ME
                disp('function evaluation failed:');
                disp(ME.message);
                stopOptimization=1;
            end
            int.mediator.setAnswer(res);
            if (createStopBox == 1) ; drawnow; end;
            if ((stopOptimization==1) && (stopOnce==1))
                disp('User interrupt requested ...');
                stopOptimize(int);
                stopOnce=0;
            end
        end
    end
    clear global JEMediator;
catch ME
    disp('Error in evaluate!');
    disp(ME.message);
    %int.mediator.quit; % just in case
    %int.mediator='';
    
    % why should this be done more than once in the end?
    %if (nontoolboxopt == 1)
    %    if (ishandle(int.boxHandle)) , close(int.boxHandle); int.boxHandle=''; end
    %    clear global stopOptimization
    %end
end
% disp('writing back results');
% write back results
int=setResultJE(int, int.mediator.getSolution());
int=setResultArrayJE(int, int.mediator.getSolutionSet());

int.mediator.quit; % just in case
int.mediator='';

% handle the case when the optimization has been called from a users script
% and not from the toolboxes parameter estimation function (which has an
% own stop button). we decide this by checking nontoolboxopt
if nontoolboxopt == 1,
    if createStopBox==1
        if (ishandle(boxHandle)) , close(boxHandle); end
    end
    clear boxHandle
    clear global stopOptimization
end
% disp('runEvalLoop done');
