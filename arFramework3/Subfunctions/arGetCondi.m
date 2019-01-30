% val = arGetCondi( m, d, condition_variable )
% Get the value of a condition in a particular dataset
%
% Usage:
%       arGetCondi( m, d, condition_variable )
%  
% Where:
%       m                       - model ID
%       d                       - data ID
%       condition_variable      - condition variable to obtain value of
%
function val = arGetCondi( m, d, condition_variable )

    global ar;

    % Is it in the condition variable list?
    ID = find( strcmp( ar.model(m).data(d).pold, condition_variable ) );

    val = '';
    if ( isempty( ID ) )
        % No? Is it in the condition struct?
        nFields = length(ar.model(m).data(d).condition);
        for jcs = 1 : nFields
            if strcmp( ar.model(m).data(d).condition(jcs).parameter, condition_variable )
                val = ar.model(m).data(d).condition(jcs).value;
            end
        end
    else
        val = ar.model(m).data(d).fp{ID};
    end
    
    

end