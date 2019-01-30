% arPrintSteadyState( [models], [mode] )
%
% Prints information about the steady state pre-equilibrations present in
% the model.
%
%    models     refers to the model IDs in the ar structure [1:length(ar.model)]
%    mode       refers to the display mode  [0]
%               0 = full 
%               1 = only condition numbers 
%               2 = only errors
% 
% Note that mode can only be specified if the model has been specified explicitly

function arPrintSteadyState( models, mode )

    global ar;
    
    if nargin < 1
        models = 1 : length( ar.model );
    end
    if nargin < 2
        mode = 0;
    end
    
    for mno = 1 : length( models )
        m = models(mno);
        if ( isfield( ar.model(m), 'ss_condition' ) )
            model = ar.model(m);
            if ( mode ~= 2 )
                linkMatrix = zeros(numel(model.condition));

                ignoreStr = cell(1, length(model.condition));
                for ssc = 1 : length( model.ss_condition )
                    linkMatrix(model.ss_condition(ssc).src, model.ss_condition(ssc).ssLink) = 1;
                    for a = 1 : numel(model.ss_condition(ssc).ssIgnore)
                        ignoreStr{model.ss_condition(ssc).src} = sprintf( '%s %s ', ignoreStr{model.ss_condition(ssc).src}, model.x{model.ss_condition(ssc).ssIgnore(a) } );
                    end
                    
                    if ( mode == 1 )
                        fprintf( '%.3d => ', model.ss_condition(ssc).src )
                        fprintf( '%.3d ', model.ss_condition(ssc).ssLink );
                        if ( length( model.ss_condition(ssc).ssIgnore ) > 0 )
                            fprintf( 'without%s', ignoreStr{model.ss_condition(ssc).src} );
                        end
                        fprintf( '\n' );
                    end
                end

                if ( mode == 0 )
                    ml = struct;
                    for d = 1 : length( model.data )
                        ml = maxLengths(ml, model.data(d).condition);
                    end

                    conditionStr = cell(1, length(model.condition));
                    for c = 1 : length( model.condition )
                        try
                            conditionStr{c} = formatCondition( model.data(model.condition(c).dLink(1)).condition, ml );
                        catch
                            conditionStr{c} = 'model default';
                        end
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
    
    if ( length( err1 ) > 0 ) %#ok
        error( '<<< FATAL ERROR >>>\nConditions %d require equilibration, but serve as equilibration for another steady state.\nThis WILL result in undefined behaviour!\nDid you forget to call arClearEvents first?\n', err1 );
    end
    if ( length( err2 ) > 0 ) %#ok
        error( '<<< FATAL ERROR >>>\nConditions %d are used twice as target in steady state equilibration.\nThis WILL result in undefined behaviour!\nDid you forget to call arClearEvents first?\n', err2 );
    end    
end