% Check for negative concentrations

function arCheckNegConc()

    global ar;
    
    for jm = 1 : length( ar.model )
        for jc = 1 : length( ar.model(jm).condition )
            found = checkNeg( jm, jc, 'x' );
            found = found || checkNeg( jm, jc, 'z' );
        end
    end
    
    if ( found == 0 )
        fprintf( 'No simulations with values smaller than -ar.config.atol were found.\n' );
    end
end

function found = checkNeg( jm, jc, field )
    global ar;

	neg = find( sum(ar.model(jm).condition(jc).([field 'ExpSimu'])<-ar.config.atol) + sum(ar.model(jm).condition(jc).([field 'FineSimu'])<-ar.config.atol) );
	dataIDs = ar.model(jm).condition(jc).dLink;
    
    % Determine which data is involved and for which state/derived variable
    involvedData = cell( 1, length( ar.model(jm).(field) ) );
    for jd = 1 : length( dataIDs )
        for jobs = 1 : length( ar.model(jm).data(dataIDs(jd)).fy )
            tokens = strsplit( ar.model(jm).data(dataIDs(jd)).fy{jobs}, {'*','/','+','-','^','(',')', ' ' } );
            if ismember(tokens, ar.model(jm).(field))
                involvedStates = find(ismember(ar.model(jm).(field), tokens));
                involvedData{involvedStates} = union( involvedData{involvedStates}, ar.model(jm).data(dataIDs(jd)).name );
            end
        end
    end

    found = 0;
	if ( ~isempty(neg) )
        found = 1;
        for jx = 1 : length( neg )
            S = unique(involvedData{jx});
            fprintf( 'm(%d).c(%d).%s(%s)', jm, jc, field, ar.model(jm).(field){jx} );
            if ~isempty(S)
                fprintf( ' involved in: \n' );
                fprintf( '  %s\n', S{:} );
            else
                fprintf( '\n' );
            end
        end
	end
end