% D = arAddDataDefInfo(RawData, DefFile, CreateData, mixingCond, opts)
%
% This function combines the data structures created by user and read from file
%
%   D                struct      data structure from .def file
%   dataD            struct      data file read from file
%   DC               struct      Data struct explicitly added by user
%   mixingCond       double      conditions show how DC is going to be used [x,x,x]

function D = arAddDataDefInfo(dataD, D, DC, mixingCond, opts)

q=0;
if ~isempty(D)
    q=1;
    arFprintf( 1, 'Combining data definitions...' );
end

if ~exist('DC','var') || ~isstruct(DC)
    DC={};
    isDC = 0;
    isDef = 0;
    isData = 0;
else
    if mixingCond(1)==0
        DC = {};
        isDC = 0;
        isDef = 0;
        isData = 0;
    else
        if  isempty(fields(DC))
            DC = {};
            isDC = 0;
            isDef = 0;
            isData = 0;
        else
            isDC=1;
            isDef = mixingCond(2);
            isData = mixingCond(3);
        end
    end
end

D = combineDataStructs(dataD, D, DC, isDC, isDef, isData);
D = orderfields(D);

if q==1
    arFprintf( 1, '[Ok]\n' );
end

if isDC    
    writeToFile(D,opts);  
end

end



% FUNCTIONS ---------------------------------------------------------------


function D = combineDataStructs(dataD, D, DC, isDC, isDef, isData)

D.extension = dataD.extension;

if isDC
    
    dataDC.path = DC.path;
    dataDC.xlstimes = DC.xlstimes;
    dataDC.xlstimevar = DC.xlstimevar;
    dataDC.xlsheader = DC.xlsheader;
    dataDC.xlsdata = DC.xlsdata;
    dataDC.xlsdataCell = DC.xlsdataCell;    
    DC = rmfield(DC,'xlstimes');
    DC = rmfield(DC,'xlstimevar');
    DC = rmfield(DC,'xlsheader');
    DC = rmfield(DC,'xlsdata');
    DC = rmfield(DC,'xlsdataCell');  
    
    % only replaces data def, keeps data has been read from data file
    if isDef==1
        D=DC;
    end
    
    % updates data def, keeping data fields (D.xls*) unchanged
    if isDef==0
        
        if ~isempty(DC.description)
            D.description = cat(1,{'* Updated by arCreatDataStruct'},D.description);
            D.description = cat(1,D.description,DC.description);
        end
        
        if ~isempty(DC.name)
            D.name = DC.name;
        end
        
        if ~isempty(DC.condition)
            D.condition = DC.condition;
        end
        
        if ~isempty(DC.doseresponse)
            D.doseresponse = DC.doseresponse;
        end
        
        if ~isempty(DC.response_parameter)
            D.response_parameter = DC.response_parameter;
        end
        
        
        %PREDICTOR / TIME
        if ~isempty(DC.t)
            D.t = DC.t;
        end
        if ~isempty(DC.tUnits)
            D.tUnits = DC.tUnits;
        end
        if ~isempty(DC.tLim)
            D.tLim = DC.tLim;
        end
        if ~isempty(DC.tLimExp)
            if ~isnan(DC.tLimExp(1))
                D.tLimExp(1)=DC.tLimExp(1);
            end
            if ~isnan(DC.tLimExp(2))
                D.tLimExp(2)=DC.tLimExp(2);
            end
            
        end
        
        
        %INPUTS
        if ~isempty(DC.u)
            
            if ~isrow(DC.u)
                DC.u = transpose(DC.u);
            end
            
            [qc,q] = ismember(DC.u,D.u);
            idxq = nonzeros(q);
            
            % update existing parameters
            D.fu(idxq) = DC.fu(qc);
            D.uNames(idxq) = DC.uNames(qc);
            
            % add new parameters
            D.u = cat(2, D.u, DC.u(qc==0));
            D.fu = cat(1, D.fu, DC.fu(qc==0));
            D.uNames = cat(2, D.uNames, DC.uNames(qc==0));
            
        end
        
        
        % OBSERVABLES & ERRORS
        if ~isempty(DC.y)
            
            if ~isrow(DC.y)
                DC.y = transpose(DC.y);
            end
            if isrow(DC.fystd)
                DC.fystd = transpose(DC.fystd);
            end
            
            [qc,q] = ismember(DC.y,D.y);
            idxq = nonzeros(q);
            
            % update existing parameters
            D.yNames(idxq) = DC.yNames(qc);
            D.yUnits(idxq,:) = DC.yUnits(qc,:);
            D.fy(idxq) = DC.fy(qc);
            D.logfitting(idxq) = DC.logfitting(qc);
            D.logplotting(idxq) = DC.logplotting(qc);
            D.normalize(idxq) = DC.normalize(qc);
            % errors
            D.fystd(idxq) = DC.fystd(qc);
            
            
            % add new parameters
            D.y = cat(2, D.y, DC.y(qc==0));
            D.yUnits = cat(1, D.yUnits, DC.yUnits(qc==0,:));
            D.yNames = cat(2, D.yNames, DC.yNames(qc==0));
            D.fy = cat(1, D.fy, DC.fy(qc==0));
            D.logfitting = cat(2, D.logfitting, DC.logfitting(qc==0));
            D.logplotting = cat(2, D.logplotting, DC.logplotting(qc==0));
            D.normalize = cat(2, D.normalize, DC.normalize(qc==0));
            %errors
            D.fystd = cat(1, D.fystd, DC.fystd(qc==0));
            
        end
        
        
        % CONDITIONS
        if ~isempty(DC.pcond)
            
            if isrow(DC.pcond)
                DC.pcond = transpose(DC.pcond);
            end
            
            [qc,q] = ismember(DC.pcond,D.pcond);
            idxq = nonzeros(q);
            
            if isfield(DC,'fpcond') && length(DC.fpcond)==length(DC.pcond)
                % update existing parameters
                D.fpcond(idxq) = DC.fpcond(qc);
                
                % add new parameters
                D.pcond = cat(1, D.pcond, DC.pcond(qc==0));
                D.fpcond = cat(1, D.fpcond, DC.fpcond(qc==0));
            else
                warning('Inconsistency in CONDITIONS created by arCreatDataStruct.')
                warning('User added conditions ignored');
            end
            
        end
        
        
        % RANDOM
        if ~isempty(DC.prand)           
            [qc,q] = ismember(DC.prand,D.prand);
            idxq = nonzeros(q);           
            % update existing parameters
            D.rand_type(idxq) = DC.rand_type(qc);            
            % add new parameters
            D.prand = cat(2, D.prand, DC.prand(qc==0));
            D.rand_type = cat(2, D.rand_type, DC.rand_type(qc==0));            
        end        
        % PARAMETERS
        % not defined yet!
    end
 
    
    
    if isData == 0        
        D.extension = dataD.extension;
        D.path = dataD.path;
        D.xlstimevar = dataD.xlstimevar;
        D.xlstimes = dataD.xlstimes;
        D.xlsheader = dataD.xlsheader;
        D.xlsdata = dataD.xlsdata;
        D.xlsdataCell = dataD.xlsdataCell;
        if isfield(dataDC,'xlstimes')
            if D.xlstimes~=dataDC.xlstimes
                error('tExp of two data sets are not match!')
            else               
                [qc,q] = ismember(dataDC.xlsheader,D.xlsheader);
                idxq = nonzeros(q);               
                % update existing data
                D.xlsdata(:,idxq) = dataDC.xlsdata(:,qc);
                D.xlsdataCell(:,idxq) = dataDC.xlsdataCell(:,qc);               
                % add new data
                D.xlsheader = cat(2, D.xlsheader, dataDC.xlsheader(qc==0));
                D.xlsdata = cat(2, D.xlsdata, dataDC.xlsdata(:,qc==0));
                D.xlsdataCell = cat(2, D.xlsdataCell, dataDC.xlsdataCell(:,qc==0));
            end            
        end
    else
        D.extension = 'xls';
        D.path = dataDC.path;
        D.xlstimes = dataDC.xlstimes;
        D.xlstimevar = dataDC.xlstimevar;
        D.xlsheader = dataDC.xlsheader;
        D.xlsdata = dataDC.xlsdata;
        D.xlsdataCell = dataDC.xlsdataCell;        
    end    
    % remove data with NaN times
    qtimesnonnan = ~isnan(D.xlstimes);
    D.xlstimes = D.xlstimes(qtimesnonnan);
    D.xlsdata = D.xlsdata(qtimesnonnan,:);
    D.xlsdataCell = D.xlsdataCell(qtimesnonnan,:);
    D.xlstimevar = D.xlstimevar;
    D.xlsheader = D.xlsheader;   
else   
    % remove data with NaN times
    qtimesnonnan = ~isnan(dataD.xlstimes);
    D.xlstimes = dataD.xlstimes(qtimesnonnan);
    D.xlsdata = dataD.xlsdata(qtimesnonnan,:);
    D.xlsdataCell = dataD.xlsdataCell(qtimesnonnan,:);
    D.xlstimevar = dataD.xlstimevar;
    D.xlsheader = dataD.xlsheader;  
end

end

function writeToFile(D,opts)
% This function creates new data struct or updated one to .def and .xls files.
% Names of the new files ende with _AutoGen. The user may use these file to
% load data independently.

arFprintf( 1, 'Generating %s.def and %s.xls files...', D.name, D.name);

if ~exist('opts','var') || isempty(opts.datapath_args)
    DataPath = 'Data/';
else
    DataPath = opts.datapath_args;
    if DataPath(end)~='/' && DataPath(end)~='\'
        DataPath = [DataPath,'/'];
    end
end

path = [D.path D.name];

fid = fopen([path '.def'],'wt');

% Write data definition file
fprintf(fid, 'DESCRIPTION\n');
for i=1:length(D.description)
fprintf(fid, '"%s"\n',D.description{i,:});
end

if D.doseresponse
   fprintf(fid, '\nPREDICTOR-DOSERESPONSE  \t %s\n', D.response_parameter); 
else
    fprintf(fid, '\nPREDICTOR\n');
end
fprintf(fid, '%s\t\t %s   "%s"  "%s"\t\t%d\t%d',D.t, D.tUnits{1}, D.tUnits{2}, D.tUnits{3}, D.tLim(1), D.tLim(2));
if ~exist('D.tLimExp','var') && sum(isnan(D.tLimExp))==0
    fprintf(fid, '\t\t%d\t%d\n',D.tLimExp(1), D.tLimExp(2));
else
    fprintf(fid, '\n');
end


fprintf(fid, '\nINPUTS\n');
for i=1:length(D.u)
fprintf(fid, '%s\t\t\t"%s"\n',D.u{i}, D.fu{i});
end


fprintf(fid, '\nOBSERVABLES\n');
for i=1:length(D.y)
fprintf(fid, '%s\t\t %s   "%s"  "%s"\t\t%d\t%d\t\t"%s"\t"%s"\n',D.y{i}, D.yUnits{i,1}, D.yUnits{i,2}, D.yUnits{i,3}, D.normalize(i), D.logfitting(i), D.fy{i}, D.yNames{i});
end 

fprintf(fid, '\nERRORS\n');
for i=1:length(D.y)
fprintf(fid, '%s\t\t\t"%s"\n',D.y{i}, D.fystd{i});
end 

if isfield(D,'subs')
    fprintf(fid, '\nSUBSTITUTIONS\n');
    for i=1:length(D.subs.from)
        fprintf(fid, '%s\t\t\t"%s"\n',D.subs.to{i}, D.subs.from{i});
    end
end

fprintf(fid, '\nCONDITIONS\n');
for i=1:length(D.pcond)
fprintf(fid, '%s\t\t\t"%s"\n',D.pcond{i}, D.fpcond{i});
end


if ~isempty(D.prand)
    fprintf(fid, '\nRANDOM\n');
    for i=1:length(D.prand)
        if D.rand_type(i)==1
            fprintf(fid, '%s\t\t\t"NORMAL"\n',D.prand{i});
        else
            fprintf(fid, '%s\t\t\t"INDEPENDENT"\n',D.prand{i});
        end
    end
end


if isfield(D,'par') && ~isempty(D.par.pExternLabels)
    fprintf(fid, '\nPARAMETERS\n');
    fprintf(fid, '#ar.pExternLabels   ar.pExtern    ar.qFitExtern    ar.qLog10Extern    ar.lbExtern    ar.ubExtern\n');
    for i=1:length(D.par.pExternLabels)
        fprintf(fid, '%s\t%d\t%d\t%d\t%d\t%d\n',D.par.pExternLabels{i}, D.par.pExtern(i), D.par.qFitExtern(i), D.par.qLog10Extern(i), D.par.lbExtern(i), D.par.ubExtern(i));
    end
end
fclose(fid);


% Write experimental data
data = cat(2, D.xlstimes, D.xlsdata);
headers = cat(2, D.xlstimevar, D.xlsheader);
table = cat(1, headers, num2cell(data));
writecell(table, [path '.xls'], 'Sheet', 1, 'Range', 'A1');

 arFprintf( 1, '[OK]\n\n' );
end
