function testFeature()

    fprintf( 2, 'TEST FOR DATA CULLING\n' );
    fprintf( 2, 'Testing whether data is appropriately culled ... ' );

    global ar;
    data = importdata('Data/filter.xls');

    % Input A and B are actually in the model. Input C is not.
    culls{1}    = {'time', @(t)(str2num(t)>160), 'input_C', @(C)(str2num(C)>0.5)};
    culls{2}    = {'time', @(t)(str2num(t)>160), 'input_C', @(C)(str2num(C)>0.5), 'input_B', @(B)(str2num(B)>0)};
    culls{3}    = {'nExpID', @(nExpID)(str2num(nExpID)<4), 'time', @(t)(str2num(t)>160), 'input_C', @(C)(str2num(C)>0.5)};
    
    arInit;
    cAr = ar;
    
    for a = 1 : numel( culls ) 
        cull = culls{a};
        arLoadModel('test');
        arLoadData( 'filter',   1, 'xls', true, 'RemoveConditions', cull );
        remaining = filterData( data, cull );
        maskTable = maskTable( data );
        pass(a) = sum(sum((remaining.'==maskTable)|isnan(maskTable))) == numel(maskTable); %#OK
    
        ar = cAr;
    end
    
    if ( sum(pass) == numel( culls ) )
        fprintf( 2, 'PASSED\n' );
    else
        pass
        error( 'DID NOT PASS TEST' );
    end
end

function remaining = filterData( data, cull )
    dataBlock = arrayfun(@(x)num2str(x), data.data, 'UniformOutput', false);
    remove = zeros( size( dataBlock, 1 ), 1 );
    for c = 1 : 2 : numel( cull )
        mask = ismember( data.colheaders, cull{c} );
        remove = remove | cellfun( cull{c+1}, dataBlock( :, mask ) );
    end
    remaining = dataBlock(~remove, :);
    remaining = cell2mat(cellfun( @(x)str2num(x), remaining, 'UniformOutput', false ));
end

function maskTable = maskTable( data )
    global ar;
    maskTable = [];
    for jd = 1 : numel( ar.model(1).data )
        tExp = ar.model(1).data(jd).tExp;
        mask = nan( numel( data.colheaders ), numel(tExp) );
        mask(1,:) = tExp;
        mask(2,:) = str2num( ar.model(1).data(jd).fprand );
        for jc = 1 : numel( ar.model(1).data(jd).condition )
            mask( ismember( data.colheaders, ar.model(1).data(jd).condition(jc).parameter ), : ) = str2num(ar.model(1).data(jd).condition(jc).value);
        end
        for jy = 1 : numel( ar.model(1).data(jd).fy )
            mask( ismember( data.colheaders, ar.model(1).data(jd).y ), : ) = ar.model(1).data(jd).yExp.';
        end
        maskTable = [maskTable , mask];
    end
end