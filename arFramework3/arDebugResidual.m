function arDebugResidual()

    global ar;
    
    maxDataLength = 0;
    maxObsLength = 0;
    for m = 1 : length( ar.model )
        for d = 1 : length( ar.model(m).data  )
            maxDataLength = max( [ ar.model(m).data(d).name maxDataLength ] );
            maxObsLength  = max( [ max(cellfun(@length, ar.model(m).data(d).y)) maxObsLength ] );
        end
    end
  
    a = 0;
    str = sprintf('Problematic residual detected\n');
    str = sprintf('%s\nErr   Model   %s    %s #dataset    #data no.\n', str, arExtendStr('Observable',maxObsLength), arExtendStr('Data file',maxDataLength) );
    
    for m = 1 : length( ar.model )
        for d = 1 : length( ar.model(m).data )
            for e = 1 : length( ar.model(m).data(d).y )
                usedIndices = find( ~isnan( ar.model(m).data(d).yExp(:,e) ) );

                nans = find( isnan( ar.model(m).data(d).res( usedIndices, e ) ) );
                infs = find( isinf( ar.model(m).data(d).res( usedIndices, e ) ) );

                if ( ~isempty(nans) )
                    a = a + 1;
                    ids = num2str(nans.');
                    str = sprintf('%sNaN   %.5d   %s    %s [%.3d]       %s\n', str, m, arExtendStr(ar.model(m).data(d).y{e},maxObsLength), arExtendStr(ar.model(m).data(d).name, maxDataLength), d, ids );
                end
                if ( ~isempty(infs) )
                    a = a + 1;
                    ids = num2str(infs.');
                    str = sprintf('%sInf   %.5d   %s    %s [%.3d]       %s\n', str, m, arExtendStr(ar.model(m).data(d).y{e},maxObsLength), arExtendStr(ar.model(m).data(d).name, maxDataLength), d, ids );
                end
            end
        end
    end
    if ( a > 0 )
        disp(str);
    end
    
end
