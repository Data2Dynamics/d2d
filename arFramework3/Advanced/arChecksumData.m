% This function runs over all ar.model.data and creates a checksum from all
% data-related fields (and ar.model.data.checkstr)

function checkstr = arChecksumData(arStruct)
global ar
if nargin==0
    arStruct = ar;
end

if ~isfield(arStruct,'model')
    checkstr = '';
else

    checkfields = {'yExp','yExpStd','qFit','qLog10','logfitting','normalize','tExp','checkstr'};
    
    checksum = [];
    for m=1:length(arStruct.model)
	if isfield( ar.model(m), 'data' )
	   for d=1:length(arStruct.model(m).data)
	        for i=1:length(checkfields)
	            if isfield(arStruct.model(m).data(d),checkfields{i})
	                val = arStruct.model(m).data(d).(checkfields{i});
	                checksum = arAddToCheckSum(val,checksum);
	            end
	        end
	    end
	end
        
        for c=1:length(arStruct.model(m).condition)
            if isfield(arStruct.model(m).condition(c),'tEvents')
                checksum = arAddToCheckSum(arStruct.model(m).condition(c).tEvents,checksum);
            end
        end
    end
    
    h = typecast(checksum.digest,'uint8');
    checkstr = dec2hex(h)';
    checkstr = checkstr(:)';
    
    clear checksum    
end