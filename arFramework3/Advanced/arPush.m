% arPush( [name/reset] )
%
% Push parameter set onto the stack.
%
%   name/reset  Name of the parameter / option stack.   [0]
%               Specify 1 to reset the stack.
%
% arPush pushes a parameter set onto the stack. The stack is a list of
% parameters that are kept for use later. The stack contains the values
% stored in ar.p, qFit, qLog10, lb, ub, type, mean and std. arPop takes the
% last parameter set from the stack and removes it from the stack.
% The stack is kept in a last in, first out format.
%
% Example(s):
%   arPush('try a new fit');
%   ar.qFit(113) = 0;
%   arFit;
%   arPop;      % Undoes the fit and undoes setting qFit 113 to zero
%
%   arPush('test1');
%   arPush('test2');
%   arPop();     % Brings the state back to test2
%   arPop();     % Brings the state back to test1
%
%   arPush('test1');
%   arPop('discard');   % Drops the last one from the stack but doesn't
%                       % bring the state back to test1.
%
% See also arPop

function arPush( reset )
    global ar;
    global arStack;
    
    if exist( 'reset', 'var' )
        if ( ~isnumeric( reset ) )
            if strcmpi( reset, 'discard' )
                reset = 1;
            else
                name = reset;
                reset = 0;
            end
        end
    else
        reset = 0;
    end    
    
    if ~exist( 'name', 'var' )
        name = '';
    end
    
    % Do we have a compatible stack?
    valid = true;
    if ( isempty(arStack) || ~isfield( arStack, 'np' ) || ( ~strcmp( arStack.checkstr, ar.checkstr ) ) )
        valid = false;
    else
        if ( length(ar.p) ~= arStack.np )
            valid = false;
        end
    end
    
    if ( ~valid )
%         the following message confuses new users
%         disp( 'The model(s) loaded are incompatible with the stored stack. Stack invalidated.' ); 
        newStack();
    end
    
    if ( reset )
        disp( 'Starting new stack' );
        newStack();
    end
       
    % Push parameter set onto the stack
    arStack.N             = arStack.N + 1;
    N                     = arStack.N;
    arStack.p(N,:)        = ar.p + 0;
    arStack.qFit(N,:)     = ar.qFit + 0;
    arStack.qLog10(N,:)   = ar.qLog10 + 0;
    arStack.lb(N,:)       = ar.lb + 0;
    arStack.ub(N,:)       = ar.ub + 0;
    arStack.type(N,:)     = ar.type + 0;       
    arStack.mean(N,:)     = ar.mean + 0;
    arStack.std(N,:)      = ar.std + 0;
    arStack.name{N}       = name;
    
    if isfield( ar.config, 'continuousSaving' )
        save( ar.config.continuousSaving, 'arStack' );
    end    
end

function newStack()
    global ar;
    global arStack;
    
    arStack             = struct;
    arStack.checkstr    = ar.checkstr;
    arStack.N           = 0;
    arStack.np          = length(ar.p);
end