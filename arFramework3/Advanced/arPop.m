% Pop parameter set
%
% arPop pops a parameter set off the stack and makes it current

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
