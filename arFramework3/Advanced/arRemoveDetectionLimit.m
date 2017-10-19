% Remove a detection limit for a specific model m and dataset d

function arRemoveDetectionLimit( m, d )
    global ar;
    
    ar.model(m).data(d).resfunction = [];
end
