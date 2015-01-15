function int = stopOptimize(int, varargin)
% Stop a running optimization. 
%   stopOptimize(JI [,'kill'])
%       If 'kill' is given as second argument, the mediator thread is
%       stopped, relevant if optimization was stopped using CTRL-C and the
%       mediator is a running zombie.

global JEMediator
    
%disp('in Stop!');
int.mp.stopOptimize;

if (nargin > 1) && (ischar(varargin{1}) && (strcmp(varargin{1},'kill')==1))
    if (~isempty(JEMediator))
        disp('killing mediator...');
        JEMediator.quit; % just in case
        JEMediator='';
        clear global JEMediator;
        clear global stopOptimization;
    else 
        disp('no mediator to kill');
    end
end