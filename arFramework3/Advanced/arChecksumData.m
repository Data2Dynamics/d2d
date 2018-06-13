% This function runs over all ar.model.data and creates a checksum from all
% data-related fields (and ar.model.data.checkstr)

function checkstr = arChecksumData

global ar

checkfields = {'yExp','yExpStd','qFit','qLog10','logfitting','normalize','tExp','checkstr'};

checksum = [];
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        for i=1:length(checkfields)
            if isfield(ar.model(m).data(d),checkfields{i})
                val = ar.model(m).data(d).(checkfields{i});
                checksum = arAddToCheckSum(val,checksum);
            end
        end
    end
    
    for c=1:length(ar.model(m).condition)
        if isfield(ar.model(m).condition(c),'tEvents')
            checksum = arAddToCheckSum(ar.model(m).condition(c).tEvents,checksum);
        end            
    end
end

h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

clear checksum

