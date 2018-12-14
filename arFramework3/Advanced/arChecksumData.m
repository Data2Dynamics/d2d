% checkstr = arChecksumData([arStruct],[saveEvaluatedFields])
% 
% This function runs over all ar.model.data and creates a checksum from all
% data-related fields (and ar.model.data.checkstr)
% 
%  arStruct         if instead of the global ar, the checksum should be
%                   evaluated for another struct, then it is provided as
%                   first argument
% 
%  saveEvaluatedFields  [false]
%                   if true, then a workspace is saved in folder Checksums
%                   containing the field which are evaluated for
%                   calculationg the checksum

function checkstr = arChecksumData(arStruct,saveEvaluatedFields)
global ar
if ~exist('arStruct','var') || isempty(arStruct)
    arStruct = ar;
end
if ~exist('saveEvaluatedFields','var') || isempty(saveEvaluatedFields)
    saveEvaluatedFields = false;
end


if ~isfield(arStruct,'model')
    checkstr = '';
else
    
    checkfields = {'yExp','yExpStd','qFit','qLog10','logfitting','normalize','tExp','checkstr'};
    
    checksum = [];
    if saveEvaluatedFields
        arCopy = struct;
    end
    
    for m=1:length(arStruct.model)
        if isfield( arStruct.model(m), 'data' )
            if saveEvaluatedFields
                modelCopy = struct;
            end
            for d=1:length(arStruct.model(m).data)
                if saveEvaluatedFields
                    dataCopy = struct;
                end
                for i=1:length(checkfields)
                    if isfield(arStruct.model(m).data(d),checkfields{i})
                        val = arStruct.model(m).data(d).(checkfields{i});
                        if saveEvaluatedFields
                            dataCopy.(checkfields{i}) = val;
                        end
                        checksum = arAddToCheckSum(val,checksum);
                    end
                end
                if saveEvaluatedFields
                    modelCopy.data(d) = dataCopy;
                end
            end
        end
               
        for c=1:length(arStruct.model(m).condition)
            cCopy = struct;
            if isfield(arStruct.model(m).condition(c),'tEvents')
                val = arStruct.model(m).condition(c).tEvents;
                cCopy.tEvents = val;
                checksum = arAddToCheckSum(val,checksum);
            end
            if saveEvaluatedFields
                condCopy(c) = cCopy;
            end
        end
        if length(arStruct.model(m).condition)>0 && saveEvaluatedFields
            modelCopy.condition = condCopy;
        end
        if saveEvaluatedFields && ~isempty(fieldnames(modelCopy))
            arCopy.model(m) = modelCopy;
        end
    end
    
    h = typecast(checksum.digest,'uint8');
    checkstr = dec2hex(h)';
    checkstr = checkstr(:)';

    if saveEvaluatedFields
        arSaveChecksumCopy(arCopy,'data',checkstr); 
    end
    
    clear checksum    
end

