% Clears condition constraints

function arClearConditionConstraints( )
    global ar;
    if isfield( ar, 'conditionconstraints' )
        ar = rmfield( ar, 'conditionconstraints' );
    end
    ar.sconstr = [];
end