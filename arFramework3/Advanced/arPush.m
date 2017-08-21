% Push parameter set
%
% arPush pushes a parameter set onto the stack

function arPush( reset )
    global ar;
    global arStack;
    
    if ~exist( 'reset', 'var' )
        reset = 0;
    else
        if ( ~isnumeric( reset ) && strcmpi( reset, 'discard' ) )
            reset = 1;
        elseif ( ~isnumeric( reset ) )
            name = reset;
        end
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
        disp( 'The model(s) loaded are incompatible with the stored stack. Stack invalidated.' );
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
end

function newStack()
    global ar;
    global arStack;
    
    arStack             = struct;
    arStack.checkstr    = ar.checkstr;
    arStack.N           = 0;
    arStack.np          = length(ar.p);
end