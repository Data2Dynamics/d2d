% status = ckWriteDataDef(file,data,force)
%  
%       file        filename or filename.def
% 
%       data      ar.model.data
% 
%       conv        output von arWriteModelDef_logX
% 
%       zeroVal     which negative value should be used instead of log(0) ?
%                   [-30]
%                   either string or number
% 
% Example:
% modelname = 'FullModel';
% conv = arWriteModelDef_logX(ar.model,['Models/',modelname,'_logX.def']);


function varargout = arWriteDataDef_logX(data,file,conv,zeroVal)
if ~exist('zeroVal','var') || isempty(zeroVal)
    zeroVal = '-30';
elseif isnumeric(zeroVal)
    zeroVal = num2str(zeroVal);
end

[a,b,c]=fileparts(file);
if isempty(c)
    file = [file,'.def'];
end

fid = fopen(file,'w');

fprintf(fid,'%s\n','DESCRIPTION');
if(isfield(data,'description'))
    for i=1:length(data.description)
        fprintf(fid,'"%s"\n',data.description{i});
    end
end
fprintf(fid,'"arWriteDataDef_LogX.m: Automatic conversion to logX %s"\n',datestr(now));
fprintf(fid,'\n');

if(data.doseresponse)   
    fprintf(fid,'%s\n','PREDICTOR-DOSERESPONSE');
    fprintf(fid,'%s\t',data.response_parameter);
else
    fprintf(fid,'%s\n','PREDICTOR');
end
if(isfield(data,'predictor'))
%     for i=1:length(data.predictor)
%         fprintf(fid,'%s\n',data.predictor{i});
%     end
else
    fprintf(fid,'%s\t%s\t%s\t%s\t%d\t%d\n',data.t,data.tUnits{1},data.tUnits{2},data.tUnits{3},data.tLim(1),data.tLim(2));    
end

fprintf(fid,'\n');
    
fprintf(fid,'%s\n','INPUTS');
if(isfield(data,'inputs'))
    for i=1:length(data.fu)
        fprintf(fid,'%s\n',data.fu{i});
    end
end
fprintf(fid,'\n');
    
fprintf(fid,'%s\n','OBSERVABLES');
for i=1:length(data.yNames)
    fprintf(fid,'%s\t%s\t%s\t%s\t%i\t%i\t"%s"\t"%s"\n',data.y{i},...
        data.yUnits{i,1},data.yUnits{i,2},data.yUnits{i,3},data.normalize(i),data.logfitting(i),data.fy{i},data.yNames{i});
end
fprintf(fid,'\n');
    

fprintf(fid,'%s\n','ERRORS');
for i=1:length(data.fystd)
    fprintf(fid,'%s\t"%s"\n',data.y{i},data.fystd{i});
end
fprintf(fid,'\n');
        

fprintf(fid,'%s\n','CONDITIONS');
% if(isfield(data,'conditions'))
%     for i=1:length(data.conditions)
%         if(ischar(data.conditions.value{i}))
%             fprintf(fid,'%s\t"%s"\n',data.conditions.name{i},data.conditions.value{i});
%         else
%             fprintf(fid,'%s\t"%f"\n',data.conditions.name{i},data.conditions.value{i});
%         end
%     end
% end

% for i=1:length(data.pold)
%     fprintf(fid,'%s\t"%s"\n',data.pold{i},data.fp{i});
% end



for i=1:length(data.fp)
    if(strcmp(data.pold{i},data.fp{i})~=1)
        
        if ~isempty(strmatch(data.pold{i},conv.statesInit,'exact'))
            if strcmp(data.fp{i},'(0)')==1
                fprintf(fid,'%s\t"%s"\n',[data.pold{i},'LG'],[zeroVal]);
            else
                fprintf(fid,'%s\t"%s"\n',[data.pold{i},'LG'],...
                    ['log(',char(arSubs(arSubs(arSym(data.fp{i}),conv.States,conv.StatesLGexp),conv.StatesInit,conv.StatesInitLGexp)),')']);
            end
        else
            fprintf(fid,'%s\t"%s"\n',...
                data.pold{i},...
                char(arSubs(arSubs(arSym(data.fp{i}),conv.States,conv.StatesLGexp),conv.StatesInit,conv.StatesInitLGexp)));
        end
        
%         fprintf(fid,'%s\t"%s"\n',data.pold{i},data.fp{i});
    end
end
fprintf(fid,'\n');


fprintf(fid,'%s\n','RANDOM');
if(isfield(data,'random'))
    for i=1:length(data.random)
        fprintf(fid,'%s\n',data.random{i});
    end
end
fprintf(fid,'\n');
    

fprintf(fid,'%s\n','PARAMETERS');
if(isfield(data,'parameters'))
    for i=1:length(data.parameters)
        fprintf(fid,'%s\n',data.parameters{i});
    end
end
fprintf(fid,'\n');


status = fclose(fid);
if(nargout>1)
    varargout{1} = status;
end

