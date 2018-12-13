% arPop( [arg] )
%
% Pops parameter set from the stack and sets it as current.
%
%   arg     'silent', don't display which parameter set we switched to.
%           'discard', discard last one from the stack but don't set it.
%
% arPop pops a parameter set onto the stack. The stack is a list of
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
% See also arPush
function arPop( discard )
    global ar;
    global arStack;
    
    silent = 0;
    if ~exist( 'discard', 'var' )
        discard = 0;
    else
        if ( ~isnumeric( discard ) )
            if strcmpi( discard, 'silent' )
                silent = 1;
            end            
            if strcmpi( discard, 'discard' )
                discard = 1;
            else
                discard = 0;
            end
        end
    end
    
    % Do we have a compatible stack?
    valid = true;
    if ( isempty(arStack) || ~isfield( arStack, 'np' ) || ( ~strcmp( arStack.checkstr, ar.checkstr ) ) )
        valid = false;
    else
        if ( length(ar.p) ~= arStack.np )
            valid = false;
        end
        if ( arStack.N < 1 )
            disp( 'No more stack left to pop' );
            return;
        end
    end
    
    if ( ~valid )
        disp( 'The model(s) loaded are incompatible with the stored stack or there is no stack' );
        return;
    end
    
    % Push parameter set onto the stack
    N                   = arStack.N;
    
    if ~discard
        ar.p            = arStack.p(N,:);
        ar.qFit         = arStack.qFit(N,:);
        ar.qLog10       = arStack.qLog10(N,:);
        ar.lb           = arStack.lb(N,:);
        ar.ub           = arStack.ub(N,:);
        ar.type         = arStack.type(N,:);       
        ar.mean         = arStack.mean(N,:);
        ar.std          = arStack.std(N,:);
        
        if ~isempty( arStack.name{N} )
            if ( ~silent )
                fprintf( 'Switched to %s\n', arStack.name{N} );
            end
        end
    end

    arStack.p(N,:)      = [];
    arStack.qFit(N,:)   = [];
    arStack.qLog10(N,:) = [];
    arStack.lb(N,:)     = [];
    arStack.ub(N,:)     = [];
    arStack.type(N,:)   = [];
    arStack.mean(N,:)   = [];
    arStack.std(N,:)    = [];
    arStack.name{N}     = [];
    
    arStack.N   = arStack.N - 1;
end
