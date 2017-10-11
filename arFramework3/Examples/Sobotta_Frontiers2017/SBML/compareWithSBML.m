q = importdata('validation_with_copasi.txt');

m = 1;
c = arFindCondition('variable');
if ( sum(ismember( q.data(:,1), ar.model(m).condition(c).tFine )) ~= numel( q.data(:,1) ) )
    ar.model(m).condition(c).tExtra = q.data(:,1);
    arLink;
end

% Uniqueize the names
mo = struct;
for a = 1 : numel( ar.model.x )
    name = strrep( ar.model(m).xNames{a}, ' ', '_' );
    if ~isfield(mo, name)
        mo.(name) = 1;
    else
        mo.(name) = mo.(name)+1;
        ar.model(m).xNames{a} = sprintf( '%s_%d', ar.model(m).xNames{a}, mo.(name) );
    end
end

time = q.data(:,1);
arSimu(false,true,true);
for a = 1 : numel( ar.model(m).x )
    curName = ar.model(m).xNames{a};
    idx = ismember( q.colheaders, curName );
    copasi = q.data( :, idx );
    this = ar.model(m).condition(c).xFineSimu( ismember( ar.model(m).condition(c).tFine, time ), a );
    
    arSubplot( [], numel(ar.model(m).x), a, 'states' );
    plot( time, copasi, 'k' ); hold on;
    plot( time, this, 'r' );
    title( curName );
end
%q.data