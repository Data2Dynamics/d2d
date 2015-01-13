function arIn = arCompress(file)
% Removes some fields from ar to save memory
% 
% if a file is provided, then the compressed struct is only saved 
% i.e. it is overwritten by the uncompressed ar after saving

if(~exist('file','var'))
    file = '';
end

global ar
if(~isempty(file) || nargout>0)
    arIn = ar;
end

if(~isfield(ar.model(1).condition(1),'suFineSimu'))
    arSimu(true,true,true);
end

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        ar.model(m).condition(c).uFineSimu_dim = size(ar.model(m).condition(c).uFineSimu);
        ar.model(m).condition(c).uFineSimu = [];
        ar.model(m).condition(c).vFineSimu_dim = size(ar.model(m).condition(c).vFineSimu);
        ar.model(m).condition(c).vFineSimu = [];
        ar.model(m).condition(c).xFineSimu_dim = size(ar.model(m).condition(c).xFineSimu);
        ar.model(m).condition(c).xFineSimu = [];
        ar.model(m).condition(c).zFineSimu_dim = size(ar.model(m).condition(c).zFineSimu);
        ar.model(m).condition(c).zFineSimu = [];
        
        ar.model(m).condition(c).suFineSimu_dim = size(ar.model(m).condition(c).suFineSimu);
        ar.model(m).condition(c).suFineSimu = [];
        ar.model(m).condition(c).svFineSimu_dim = size(ar.model(m).condition(c).svFineSimu);
        ar.model(m).condition(c).svFineSimu = [];
        ar.model(m).condition(c).sxFineSimu_dim = size(ar.model(m).condition(c).sxFineSimu);
        ar.model(m).condition(c).sxFineSimu = [];
        ar.model(m).condition(c).szFineSimu_dim = size(ar.model(m).condition(c).szFineSimu);
        ar.model(m).condition(c).szFineSimu = [];
    end
end

ar.isCompressed = 1;

if(~isempty(file))
    if(~isempty(ar.config.savepath))
        pfad = ar.config.savepath;
    else
        pfad = '.';
    end
    [pfad,'/',file]
    save([pfad,'/',file],'ar');
    ar = arIn;
end


