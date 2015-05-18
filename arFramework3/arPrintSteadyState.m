% arPrintSteadyState
%
% Prints information about the steady state pre-equilibrations present in
% the model.
%
% Usage:
%    arPrintSteadyState( (model) )
%
% Optional argument:
%    Model refers to the model number in the ar structure. If omitted, all
%    models are traversed.

function arPrintSteadyState( models )

    global ar;
    
    if nargin < 1
        models = 1 : length( ar.model );
    end
    
    for mno = 1 : length( models )
        m = models(mno);
        if ( isfield( ar.model(m), 'ss_condition' ) )
            model = ar.model(m);
            linkMatrix = zeros(numel(model.condition));

            ignoreStr = cell(1, length(model.ss_condition));
            for ssc = 1 : length( model.ss_condition )
                linkMatrix(model.ss_condition(ssc).src, model.ss_condition(ssc).ssLink) = 1;
                for a = 1 : numel(model.ss_condition(ssc).src)
                    ignoreStr{model.ss_condition(ssc).src(a)} = sprintf( '%s ', model.x{model.ss_condition(ssc).ssIgnore } );
                end
            end

            ml = struct;
            for d = 1 : length( model.data )
                ml = maxLengths(ml, model.data(d).condition);
            end

            conditionStr = cell(1, length(model.condition));
            for c = 1 : length( model.condition )
                conditionStr{c} = formatCondition( model.data(model.condition(c).dLink(1)).condition, ml );
            end

            maxStrLen = max(cellfun(@length, conditionStr));
            fprintf( '#m  | #c  | %s | #c  | %s | Non-equilibrated states\n', arExtendStr( 'Target condition variables', maxStrLen ), arExtendStr( 'Steady state condition variables', maxStrLen ) );

            for to = 1 : length( model.condition )
                printed = 0;
                for from = 1 : length( model.condition )
                    if ( linkMatrix( from, to ) )
                        fprintf( '%.3d | %.3d | %s | %.3d | %s | %s \n', m, to, arExtendStr( conditionStr{to}, maxStrLen ), from, arExtendStr( conditionStr{from}, maxStrLen ), ignoreStr{from} );
                        printed = 1;
                    end
                end
                % No steady state equilibration linked up
                if ( ~printed )
                    fprintf( '%.3d | %.3d | %s | No equilibration assigned\n', m, to, arExtendStr( conditionStr{to}, maxStrLen ) );
                end
            end

            sanityCheck(model);
        else
            fprintf( 'No steady states specified.\n' );
        end
    end
end

function ml = maxLengths( ml, condPars )
    for a = 1 : length( condPars )
        if ( ~isfield( ml, condPars(a).parameter ) )
            ml.(condPars(a).parameter) = numel(condPars(a).value);
        else
            ml.(condPars(a).parameter) = max( ml.(condPars(a).parameter), numel(condPars(a).value) );
        end
    end
end

function str = formatCondition( condPars, ml )
    str = [];
    
    p = {condPars.parameter};
    [~,I] = sort(cellfun(@length,p));

    for a = 1 : length( condPars )
        str = sprintf( '%s%s = %s : ', str, condPars(I(a)).parameter, arExtendStr( condPars(I(a)).value, ml.(condPars(I(a)).parameter) ) );
    end
    
    str = str(1:end-2);
end

function sanityCheck( model )
    ss_sources = zeros(numel(model.condition),1);
    ss_targets = zeros(numel(model.condition),1);
    for ssc = 1 : length( model.ss_condition )
        ss_sources( model.ss_condition(ssc).src ) = ss_sources( model.ss_condition(ssc).src ) + 1;
        ss_targets( model.ss_condition(ssc).ssLink ) = ss_targets( model.ss_condition(ssc).ssLink ) + 1;
    end
    
    % Conditions that are used both as source as well as target
    err1 = find( ( ss_sources > 0 ) & ( ss_targets > 0 ) & ( ss_sources ~= ss_targets ) );
    
    % Conditions that are targetted twice
    err2 = find( ss_targets > 1 );
    
    if ( length( err1 ) > 0 )
        fprintf( '\n<<< WARNING >>>\nConditions ' );
        fprintf( '%d ', err1 );
        fprintf( 'require equilibration, but serve as equilibration for another steady state.\nThis WILL result in undefined behaviour!\n' );
    end
    if ( length( err2 ) > 0 )
        fprintf( '\n<<< WARNING >>>\nConditions ' );
        fprintf( '%d ', err2 );
        fprintf( 'are used twice as target in steady state equilibration.\nThis WILL result in undefined behaviour!\n' );
    end    
end