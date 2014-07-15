% Load data set to next free slot
%
% arLoadData(name, m, extension, removeEmptyObs)
%
% name              filename of data definition file
% m                 target position (int) for model                [last loaded model]
%                   or: model name (string)
% extension         data file name-extension: 'xls', 'csv'          ['xls']
%                   'none' = don't load data                           
% removeEmptyObs    remove observation without data                [false]
%
%
% In der ersten Spalte:
% 1)    Die Messzeitpunkte (duerfen mehrfach vorkommen).
%
% Danach in Spalten beliebiger Reihenfolge:
% 2)    Die experimentellen Bedingungen (z.B. "input_IL6" und "input_IL1").
% 3)    Die Datenpunkte fuer die einzelnen Spezies (z.B. "P_p38_rel").
%
% Anmerkungen:
% 1)    In Spaltenkoepfen duerfen keine Zeichen vorkommen die mathematischen
%       Operationen entsprechen (z.B. "-" oder "+").
% 2)    Ich habe Stimulationen immer den Praefix "input_" gegeben. Bei Spezies
%       bedeuten die Suffixe "_rel" und "_au": relative Phosphorylierung und
%       arbitrary units, je nachdem.
%
% Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)

function arLoadData(name, m, extension, removeEmptyObs, dpPerShoot)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

% load model from mat-file
if(~exist('Data','dir'))
    error('folder Data/ does not exist')
end
if(~exist(['Data/' name '.def'],'file'))
    error('data definition file %s.def does not exist in folder Data/', name)
end

if(~exist('m','var') || isempty(m))
    m = length(ar.model);
end
if(exist('m','var') && ischar(m))
    for jm=1:length(ar.model)
        if(strcmp(m, ar.model(jm).name))
            m = jm;
        end
    end
    if(ischar(m))
        error('Model %s was not found', m);
    end
end

if(exist('extension','var') && isnumeric(extension) && ...
        ~(isempty(extension) && nargin>3))
    error(['arLoadData(name, m, d, ...) input argument d is deprecated !!! ' ...
        'Please see new usage arLoadModel(name, m, extension, removeEmptyObs) and function help text.']);
end

if(isfield(ar.model(m), 'data'))
    d = length(ar.model(m).data) + 1;
else
    ar.model(m).data = [];
    d = 1;
end

if(~exist('extension','var') || isempty(extension))
    extension = 'xls';
end
if(~exist('removeEmptyObs','var'))
    removeEmptyObs = false;
else
    if(ischar(removeEmptyObs))
        error(['arLoadData(name, m, d, ...) input argument d is deprecated !!! ' ...
            'Please see new usage arLoadModel(name, m, extension, removeEmptyObs) and function help text.']);
    end
end

if(exist('dpPerShoot','var') && dpPerShoot>0)
    if(~isfield(ar,'ms_count_snips'))
        ar.model(m).ms_count = 0;
        ar.ms_count_snips = 0;
        ar.ms_strength = 0;
        ar.ms_threshold = 1e-5;
        ar.ms_violation = [];
    end
else
    dpPerShoot = 0;
end

% initial setup
ar.model(m).data(d).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');

fprintf('\nloading data #%i, from file Data/%s.def...', d, name);
fid = fopen(['Data/' name '.def'], 'r');

% DESCRIPTION
str = textscan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(~strcmp(str{1},'DESCRIPTION'))
    error('parsing data %s for DESCRIPTION', name);
end

% check version
if(strcmp(str{1},'DESCRIPTION'))
    % def_version = 1;
elseif(strcmp(str{1},'DESCRIPTION-V2'))
    error('DESCRIPTION-V2 not supported yet');
else
    error('invalid version identifier: %s', cell2mat(str{1}));
end

% read comments
str = textscan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
ar.model(m).data(d).description = {};
while(~strcmp(str{1},'PREDICTOR') && ~strcmp(str{1},'PREDICTOR-DOSERESPONSE'))
    ar.model(m).data(d).description(end+1,1) = str{1}; %#ok<*AGROW>
    str = textscan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
end

% PREDICTOR
if(strcmp(str{1},'PREDICTOR-DOSERESPONSE'))
    ar.model(m).data(d).doseresponse = true;
    str = textscan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
    ar.model(m).data(d).response_parameter = cell2mat(str{1});
    fprintf('dose-response to %s\n', ar.model(m).data(d).response_parameter);
else
    ar.model(m).data(d).doseresponse = false;
    ar.model(m).data(d).response_parameter = '';
    fprintf('\n');
end
C = textscan(fid, '%s %s %q %q %n %n %n %n\n',1, 'CommentStyle', ar.config.comment_string);
ar.model(m).data(d).t = cell2mat(C{1});
ar.model(m).data(d).tUnits(1) = C{2};
ar.model(m).data(d).tUnits(2) = C{3};
ar.model(m).data(d).tUnits(3) = C{4};
ar.model(m).data(d).tLim = [C{5} C{6}];
ar.model(m).data(d).tLimExp = [C{7} C{8}];
if(isnan(ar.model(m).tLim(1)))
    ar.model(m).tLim(1) = 0;
end
if(isnan(ar.model(m).tLim(2)))
    ar.model(m).tLim(2) = 10;
end
if(isnan(ar.model(m).data(d).tLimExp(1)))
    ar.model(m).data(d).tLimExp(1) = ar.model(m).tLim(1);
end
if(isnan(ar.model(m).data(d).tLimExp(2)))
    ar.model(m).data(d).tLimExp(2) = ar.model(m).tLim(2);
end

% INPUTS
str = textscan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(~strcmp(str{1},'INPUTS'))
    error('parsing data %s for INPUTS', name);
end
C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
ar.model(m).data(d).fu = ar.model(m).fu;
while(~strcmp(C{1},'OBSERVABLES'))
    qu = ismember(ar.model(m).u, C{1}); %R2013a compatible
    if(sum(qu)~=1)
        error('unknown input %s', cell2mat(C{1}));
    end
    ar.model(m).data(d).fu(qu) = C{2};
    
    C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
end

% input parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fu, 'UniformOutput', false);
ar.model(m).data(d).pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}); %R2013a compatible

% OBSERVABLES
if(isfield(ar.model(m),'y'))
    ar.model(m).data(d).y = ar.model(m).y;
    ar.model(m).data(d).yNames = ar.model(m).yNames;
    ar.model(m).data(d).yUnits = ar.model(m).yUnits;
    ar.model(m).data(d).normalize = ar.model(m).normalize;
    ar.model(m).data(d).logfitting = ar.model(m).logfitting;
    ar.model(m).data(d).logplotting = ar.model(m).logplotting;
    ar.model(m).data(d).fy = ar.model(m).fy;
else 
    ar.model(m).data(d).y = {};
    ar.model(m).data(d).yNames = {};
    ar.model(m).data(d).yUnits = {};
    ar.model(m).data(d).normalize = [];
    ar.model(m).data(d).logfitting = [];
    ar.model(m).data(d).logplotting = [];
    ar.model(m).data(d).fy = {};
end
C = textscan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'ERRORS'))
    qyindex = ismember(ar.model(m).data(d).y, C{1});
    if(sum(qyindex)==1)
        yindex = find(qyindex);
    elseif(sum(qyindex)==0)
        yindex = length(ar.model(m).data(d).y) + 1;
    else
        error('multiple matches for %s', cell2mat(C{1}))
    end
    ar.model(m).data(d).y(yindex) = C{1};
    ar.model(m).data(d).yUnits(yindex,1) = C{2};
    ar.model(m).data(d).yUnits(yindex,2) = C{3};
    ar.model(m).data(d).yUnits(yindex,3) = C{4};
    ar.model(m).data(d).normalize(yindex) = C{5};
    ar.model(m).data(d).logfitting(yindex) = C{6};
    ar.model(m).data(d).logplotting(yindex) = C{6};
    ar.model(m).data(d).fy(yindex,1) = C{7};
    if(~isempty(cell2mat(C{8})))
        ar.model(m).data(d).yNames(yindex) = C{8};
    else
        ar.model(m).data(d).yNames(yindex) = ar.model(m).data(d).y(yindex);
    end
    C = textscan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', ar.config.comment_string);
    if(sum(ismember(ar.model(m).x, ar.model(m).data(d).y{yindex}))>0) %R2013a compatible
        error('%s already defined in STATES', ar.model(m).data(d).y{yindex});
    end
    if(sum(ismember(ar.model(m).u, ar.model(m).data(d).y{end}))>0) %R2013a compatible
        error('%s already defined in INPUTS', ar.model(m).data(d).y{end});
    end
    if(sum(ismember(ar.model(m).z, ar.model(m).data(d).y{end}))>0) %R2013a compatible
        error('%s already defined in DERIVED', ar.model(m).data(d).y{end});
    end
end

% observation parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fy, 'UniformOutput', false);
ar.model(m).data(d).py = setdiff(setdiff(vertcat(varlist{:}), union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
for j=1:length(ar.model(m).data(d).fy)
    varlist = symvar(ar.model(m).data(d).fy{j});
    ar.model(m).data(d).py_sep(j).pars = setdiff(setdiff(varlist, union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
    
    % exclude parameters form model definition
    ar.model(m).data(d).py_sep(j).pars = setdiff(ar.model(m).data(d).py_sep(j).pars, ar.model(m).px);
end

% ERRORS
if(isfield(ar.model(m),'y'))
    ar.model(m).data(d).fystd = ar.model(m).fystd;
else
    ar.model(m).data(d).fystd = cell(size(ar.model(m).data(d).fy));
end
C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'INVARIANTS') && ~strcmp(C{1},'DERIVED') && ~strcmp(C{1},'CONDITIONS'))
    qyindex = ismember(ar.model(m).data(d).y, C{1});
    if(sum(qyindex)==1)
        yindex = find(qyindex);
    elseif(sum(qyindex)==0)
        yindex = length(ar.model(m).data(d).y) + 1;
    else
        error('multiple matches for %s', cell2mat(C{1}))
    end
    ar.model(m).data(d).fystd(yindex) = C{2};
    C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
end

if(length(ar.model(m).data(d).fystd)<length(ar.model(m).data(d).fy))
    error('some observables do not have an error model defined');
end

% error parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fystd, 'UniformOutput', false);
ar.model(m).data(d).pystd = setdiff(vertcat(varlist{:}), union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ... %R2013a compatible
    ar.model(m).data(d).y), ar.model(m).t));
for j=1:length(ar.model(m).data(d).fystd)
    varlist = symvar(ar.model(m).data(d).fystd{j});
	ar.model(m).data(d).py_sep(j).pars = union(ar.model(m).data(d).py_sep(j).pars, ... %R2013a compatible
        setdiff(varlist, union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ar.model(m).data(d).y))); %R2013a compatible
    
    % exclude parameters form model definition
    ar.model(m).data(d).py_sep(j).pars = setdiff(ar.model(m).data(d).py_sep(j).pars, ar.model(m).px);
end

% DERIVED
if(strcmp(C{1},'DERIVED'))
    error(['There is no need for a section DERIVED in data definition file! ' ...
        'Please remove and see usage in: ' ...
        'https://bitbucket.org/d2d-development/d2d-software/wiki/Setting%20up%20models']);
end
% INVARIANTS
if(strcmp(C{1},'INVARIANTS'))
    error(['Section INVARIANTS in data definition file is deprecated! ' ...
        'Please remove and see usage in: ' ...
        'https://bitbucket.org/d2d-development/d2d-software/wiki/Setting%20up%20models']);
end

% collect parameters needed for OBS
ptmp = union(ar.model(m).px, ar.model(m).pu);
ar.model(m).data(d).p = union(ptmp, union(ar.model(m).data(d).pu, ar.model(m).data(d).py)); %R2013a compatible
ar.model(m).data(d).p = union(ar.model(m).data(d).p, ar.model(m).data(d).pystd); %R2013a compatible

% replace filename
ar.model(m).data(d).p = strrep(ar.model(m).data(d).p, '_filename', ['_' ar.model(m).data(d).name]);
ar.model(m).data(d).fy = strrep(ar.model(m).data(d).fy, '_filename', ['_' ar.model(m).data(d).name]);
ar.model(m).data(d).fystd = strrep(ar.model(m).data(d).fystd, '_filename', ['_' ar.model(m).data(d).name]);

% CONDITIONS
C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
ar.model(m).data(d).fp = transpose(ar.model(m).data(d).p);
ptmp = ar.model(m).p;
qcondparamodel = ismember(ar.model(m).data(d).p, strrep(ptmp, '_filename', ['_' ar.model(m).data(d).name])); %R2013a compatible
qmodelparacond = ismember(strrep(ptmp, '_filename', ['_' ar.model(m).data(d).name]), ar.model(m).data(d).p); %R2013a compatible
ar.model(m).data(d).fp(qcondparamodel) = strrep(ar.model(m).fp(qmodelparacond), '_filename', ['_' ar.model(m).data(d).name]);
while(~isempty(C{1}) && ~strcmp(C{1},'RANDOM'))
    qcondpara = ismember(ar.model(m).data(d).p, C{1}); %R2013a compatible
    if(sum(qcondpara)>0)
        ar.model(m).data(d).fp{qcondpara} = ['(' cell2mat(C{2}) ')'];
    else
        warning('unknown parameter in conditions %s', cell2mat(C{1}));
    end
    C = textscan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
end

% extra conditional parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fp, 'UniformOutput', false);
ar.model(m).data(d).pcond = setdiff(vertcat(varlist{:}), ar.model(m).data(d).p); %R2013a compatible
      
% collect parameters conditions
pcond = union(ar.model(m).data(d).p, ar.model(m).data(d).pcond); %R2013a compatible

% RANDOM
ar.model(m).data(d).prand = {};
ar.model(m).data(d).rand_type = [];
C = textscan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
while(~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS'))
    ar.model(m).data(d).prand{end+1} = cell2mat(C{1});
    if(strcmp(C{2}, 'INDEPENDENT'))
        ar.model(m).data(d).rand_type(end+1) = 0;
    else
        warning('unknown random type %s', cell2mat(C{2}));  %#ok<WNTAG>
    end
    C = textscan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
end

% PARAMETERS
if(~isfield(ar, 'pExternLabels'))
    ar.pExternLabels = {};
    ar.pExtern = [];
    ar.qFitExtern = [];
    ar.qLog10Extern = [];
    ar.lbExtern = [];
    ar.ubExtern = [];
end
C = textscan(fid, '%s %f %n %n %n %n\n',1, 'CommentStyle', ar.config.comment_string);
while(~isempty(C{1}))
    ar.pExternLabels(end+1) = C{1};
    ar.pExtern(end+1) = C{2};
    ar.qFitExtern(end+1) = C{3};
    ar.qLog10Extern(end+1) = C{4};
    ar.lbExtern(end+1) = C{5};
    ar.ubExtern(end+1) = C{6};
    C = textscan(fid, '%s %f %n %n %n %n\n',1, 'CommentStyle', ar.config.comment_string);
end

% plot setup
if(isfield(ar.model(m).data(d), 'response_parameter') && ...
        ~isempty(ar.model(m).data(d).response_parameter))
    if(sum(ismember(ar.model(m).data(d).p ,ar.model(m).data(d).response_parameter))==0 && ... %R2013a compatible
            sum(ismember(ar.model(m).data(d).pcond ,ar.model(m).data(d).response_parameter))==0) %R2013a compatible
        error('invalid response parameter %s', ar.model(m).data(d).response_parameter);
    end
end
if(~isfield(ar.model(m), 'plot'))
    ar.model(m).plot(1).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
else
    ar.model(m).plot(end+1).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
end
ar.model(m).plot(end).doseresponse = ar.model(m).data(d).doseresponse;
ar.model(m).plot(end).dLink = d;
ar.model(m).plot(end).ny = length(ar.model(m).data(d).y);
ar.model(m).plot(end).condition = {};
jplot = length(ar.model(m).plot);

fclose(fid);

% XLS file
if(~strcmp(extension,'none') && ((exist(['Data/' name '.xls'],'file') && strcmp(extension,'xls')) || ...
        (exist(['Data/' name '.csv'],'file') && strcmp(extension,'csv'))))
    fprintf('loading data #%i, from file Data/%s.%s...\n', d, name, extension);
    
    % read from file
    if(strcmp(extension,'xls'))
        warntmp = warning;
        warning('off','all')
        [data, Cstr] = xlsread(['Data/' name '.xls']);
        if(length(data(1,:))>length(Cstr(1,:)))
            data = data(:,1:length(Cstr(1,:)));
        end
        
        warning(warntmp);
        
        header = Cstr(1,2:end);
        times = data(:,1);
        qtimesnonnan = ~isnan(times);
        times = times(qtimesnonnan);
        data = data(qtimesnonnan,2:end);
        if(size(data,2)<length(header))
            data = [data nan(size(data,1),length(header)-size(data,2))];
        end
        
        Cstr = Cstr(2:end,2:end);
        dataCell = cell(size(data));
        for j1 = 1:size(data,1)
            for j2 = 1:size(data,2)
                if(isnan(data(j1,j2)))
                    if(j1<=size(Cstr,1) && j2<=size(Cstr,2) && ~isempty(Cstr{j1,j2}))
                        dataCell{j1,j2} = Cstr{j1,j2};
                    else
                        dataCell{j1,j2} = header{j2};
                    end
                else
                    dataCell{j1,j2} = num2str(data(j1,j2));
                end
            end
        end
        
    elseif(strcmp(extension,'csv'))
        fid = fopen(['Data/' name '.csv'], 'r');
        
        C = textscan(fid,'%s\n',1);
        C = textscan(C{1}{1},'%q','Delimiter',',');
        C = C{1};
        header = C(2:end)';
        
        data = nan(0, length(header));
        dataCell = cell(0, length(header));
        times = [];
        rcount = 1;
        C = textscan(fid,'%q',length(header)+1,'Delimiter',',');
        while(~isempty(C{1}))
            C = strrep(C{1}',',','.');
            times(rcount,1) = str2double(C{1});
            for j=1:length(header)
                data(rcount,j) = str2double(C{j+1});
                dataCell{rcount,j} = C{j+1};
            end
            C = textscan(fid,'%q',length(header)+1,'Delimiter',',');
            rcount = rcount + 1;
        end
        
        fclose(fid);
    end
    
    % random effects
    qrandis = ismember(header, ar.model(m).data(d).prand); %R2013a compatible
    if(sum(qrandis) > 0)
        qobs = ismember(header, ar.model(m).data(d).y); %R2013a compatible
        
        randis_header = header(qrandis);
        if ~isempty(dataCell)
            [randis, ~, jrandis] = uniqueRowsCA(dataCell(:,qrandis));
        else
            [randis, ~, jrandis] = unique(data(:,qrandis),'rows');
            randis = cellstr(num2str(randis));
        end
               
        for j=1:size(randis,1)
            qvals = jrandis == j;
            tmpdata = data(qvals,qobs);
            if(sum(~isnan(tmpdata(:)))>0)
                fprintf('local random effect #%i:\n', j)
                
                if(j < size(randis,1))
                    ar.model(m).data(d+1) = ar.model(m).data(d);
                    ar.model(m).plot(jplot+1) = ar.model(m).plot(jplot);
                end
                
                for jj=1:size(randis,2)
                    fprintf('\t%20s = %s\n', randis_header{jj}, randis{j,jj})
                    
                    ar.model(m).plot(jplot).name = [ar.model(m).plot(jplot).name '_' ...
                        randis_header{jj} randis{j,jj}];
                    
                    ar.model(m).data(d).name = [ar.model(m).data(d).name '_' ...
                        randis_header{jj} randis{j,jj}];
                    
                    ar.model(m).data(d).fy = strrep(ar.model(m).data(d).fy, ...
                        randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    ar.model(m).data(d).py = strrep(ar.model(m).data(d).py, ...
                        randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    
                    ar.model(m).data(d).fystd = strrep(ar.model(m).data(d).fystd, ...
                        randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    ar.model(m).data(d).pystd = strrep(ar.model(m).data(d).pystd, ...
                        randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    
%                     ar.model(m).data(d).p = strrep(ar.model(m).data(d).p, ...
%                         randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    ar.model(m).data(d).fp = strrep(ar.model(m).data(d).fp, ...
                        randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    ar.model(m).data(d).pcond = strrep(ar.model(m).data(d).pcond, ...
                        randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    
                    for jjj=1:length(ar.model(m).data(d).py_sep)
                        ar.model(m).data(d).py_sep(jjj).pars = strrep(ar.model(m).data(d).py_sep(jjj).pars, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    end
                end
                
                if ~isempty(dataCell)
                    [ar,d] = setConditions(ar, m, d, jplot, header, times(qvals), data(qvals,:), dataCell(qvals,:), ...
                        strrep(pcond, randis_header{jj}, [randis_header{jj} randis{j,jj}]), removeEmptyObs, dpPerShoot);
                else
                    [ar,d] = setConditions(ar, m, d, jplot, header, times(qvals), data(qvals,:), dataCell, ...
                        strrep(pcond, randis_header{jj}, [randis_header{jj} randis{j,jj}]), removeEmptyObs, dpPerShoot);
                end
                if(j < size(randis,1))
                    d = d + 1;
                    jplot = jplot + 1;
                    ar.model(m).plot(jplot).dLink = d;
                end
            else
                fprintf('local random effect #%i: no matching data, skipped\n', j);
            end
        end
    else
        ar = setConditions(ar, m, d, jplot, header, times, data, dataCell, pcond, removeEmptyObs, dpPerShoot);
    end
else
    ar.model(m).data(d).condition = [];
end



function [ar,d] = setConditions(ar, m, d, jplot, header, times, data, dataCell, pcond, removeEmptyObs, dpPerShoot)

% matVer = ver('MATLAB');

% normalization of columns
nfactor = max(data, [], 1);

qobs = ismember(header, ar.model(m).data(d).y) & sum(~isnan(data),1)>0; %R2013a compatible
qhasdata = ismember(ar.model(m).data(d).y, header(qobs)); %R2013a compatible

% conditions
qcond = ismember(header, pcond); %R2013a compatible
if(sum(qcond) > 0)
    condi_header = header(qcond);
    if ~isempty(dataCell)
        [condis, ~, jcondis] = uniqueRowsCA(dataCell(:,qcond));
    else
        [condis, ~, jcondis] = unique(data(:,qcond),'rows');
        condis = mymat2cell(condis);
    end
    
    active_condi = false(size(condis(1,:)));
    tmpcondi = condis(1,:);
    for j1=2:size(condis,1)
        for j2=1:size(condis,2)
            active_condi(j2) = active_condi(j2) | (~strcmp(tmpcondi{j2}, condis{j1,j2}));
        end
    end
    
    for j=1:size(condis,1)
        
        fprintf('local condition #%i:\n', j)
        
        if(j < size(condis,1))
            if(length(ar.model(m).data) > d)
                ar.model(m).data(d+2) = ar.model(m).data(d+1);
            end
            ar.model(m).data(d+1) = ar.model(m).data(d);
        end
        
        % remove obs without data
        if(removeEmptyObs)
            for jj=find(~qhasdata)
                fprintf('\t%20s no data, removed\n', ar.model(m).data(d).y{jj});
                jjjs = find(ismember(ar.model(m).data(d).p, ar.model(m).data(d).py_sep(jj).pars)); %R2013a compatible
                jjjs = jjjs(:)';
                for jjj=jjjs
                    remove = 1;
                    for jjjj = find(qhasdata)
                        if sum(ismember(ar.model(m).data(d).py_sep(jjjj).pars, ar.model(m).data(d).p(jjj))) > 0 %R2013a compatible
                            remove = 0;
                        end
                    end
                    if remove
                        ar.model(m).data(d).fp{jjj} = '0';
                    end
                end
            end
            ar.model(m).data(d).y = ar.model(m).data(d).y(qhasdata);
            ar.model(m).data(d).yNames = ar.model(m).data(d).yNames(qhasdata);
            ar.model(m).data(d).yUnits = ar.model(m).data(d).yUnits(qhasdata,:);
            ar.model(m).data(d).normalize = ar.model(m).data(d).normalize(qhasdata);
            ar.model(m).data(d).logfitting = ar.model(m).data(d).logfitting(qhasdata);
            ar.model(m).data(d).logplotting = ar.model(m).data(d).logplotting(qhasdata);
            ar.model(m).data(d).fy = ar.model(m).data(d).fy(qhasdata);
            ar.model(m).data(d).fystd = ar.model(m).data(d).fystd(qhasdata);
            ar.model(m).data(d).py_sep = ar.model(m).data(d).py_sep(qhasdata);
        end
        
        for jj=1:size(condis,2)
            if(~isempty(condis{j,jj}))
                fprintf('\t%20s = %s\n', condi_header{jj}, condis{j,jj})
                
                qcondjj = ismember(ar.model(m).data(d).p, condi_header{jj}); %R2013a compatible
                if(sum(qcondjj)>0)
                    ar.model(m).data(d).fp{qcondjj} =  ['(' condis{j,jj} ')'];
                end
                qcondjj = ~strcmp(ar.model(m).data(d).p, ar.model(m).data(d).fp');
                if(~isnan(str2double(condis{j,jj})))
                    ar.model(m).data(d).fp(qcondjj) = strrep(ar.model(m).data(d).fp(qcondjj), ...
                        condi_header{jj}, condis{j,jj});
                end
                
                ar.model(m).data(d).condition(jj).parameter = condi_header{jj};
                ar.model(m).data(d).condition(jj).value = condis{j,jj};
                
                % plot
                if(active_condi(jj))
                    if(ar.model(m).data(d).doseresponse==0 || ~strcmp(condi_header{jj}, ar.model(m).data(d).response_parameter))
                        if(length(ar.model(m).plot(jplot).condition) >= j && ~isempty(ar.model(m).plot(jplot).condition{j}))
                            ar.model(m).plot(jplot).condition{j} = [ar.model(m).plot(jplot).condition{j} ' & ' ...
                                ar.model(m).data(d).condition(jj).parameter '=' ...
                                ar.model(m).data(d).condition(jj).value];
                        else
                            ar.model(m).plot(jplot).condition{j} = [ar.model(m).data(d).condition(jj).parameter '=' ...
                                ar.model(m).data(d).condition(jj).value];
                        end
                    end
                end
            end
        end
        
        qvals = jcondis == j;
        ar = setValues(ar, m, d, header, nfactor, data(qvals,:), times(qvals));
        ar.model(m).data(d).tLim(2) = round(max(times)*1.1);
        
        if(dpPerShoot~=0)
            [ar,d] = doMS(ar,m,d,jplot,dpPerShoot);
        end
        
        if(j < size(condis,1))
            d = d + 1;
            ar.model(m).plot(jplot).dLink(end+1) = d;
        end
    end
else
    ar.model(m).data(d).condition = [];
    
    % remove obs without data
    if(removeEmptyObs)
        for jj=find(~qhasdata)
            fprintf('\t%20s no data, removed\n', ar.model(m).data(d).y{jj});
            jjjs = find(ismember(ar.model(m).data(d).p, ar.model(m).data(d).py_sep(jj).pars)); %R2013a compatible
            jjjs = jjjs(:)';
            for jjj=jjjs
                remove = 1;
                for jjjj = find(qhasdata)
                    if sum(ismember(ar.model(m).data(d).py_sep(jjjj).pars, ar.model(m).data(d).p(jjj))) > 0 %R2013a compatible
                        remove = 0;
                    end
                end
                if(remove==1)
                    ar.model(m).data(d).fp{jjj} = '0';
                end
            end
        end
        ar.model(m).data(d).y = ar.model(m).data(d).y(qhasdata);
        ar.model(m).data(d).yNames = ar.model(m).data(d).yNames(qhasdata);
        ar.model(m).data(d).yUnits = ar.model(m).data(d).yUnits(qhasdata,:);
        ar.model(m).data(d).normalize = ar.model(m).data(d).normalize(qhasdata);
        ar.model(m).data(d).logfitting = ar.model(m).data(d).logfitting(qhasdata);
        ar.model(m).data(d).logplotting = ar.model(m).data(d).logplotting(qhasdata);
        ar.model(m).data(d).fy = ar.model(m).data(d).fy(qhasdata);
        ar.model(m).data(d).fystd = ar.model(m).data(d).fystd(qhasdata);
        ar.model(m).data(d).py_sep = ar.model(m).data(d).py_sep(qhasdata);
    end
    
    ar = setValues(ar, m, d, header, nfactor, data, times);
    
    if(dpPerShoot~=0)
        [ar,d] = doMS(ar,m,d,jplot,dpPerShoot);
    end
end


function C = mymat2cell(D)
C = cell(size(D));
for j=1:size(D,1)
    for jj=1:size(D,2)
        C{j,jj} = num2str(D(j,jj));
    end
end

function [ar,d] = doMS(ar,m,d,jplot,dpPerShoot)

tExp = ar.model(m).data(d).tExp;

if(dpPerShoot ~= 1)
    nints = ceil(length(tExp) / dpPerShoot);
    tboarders = linspace(min(tExp),max(tExp),nints+1);
else
    tboarders = union(0,tExp); %R2013a compatible
    nints = length(tboarders)-1;
end

if(nints==1)
    return;
end

fprintf('using %i shooting intervals\n', nints);
ar.model(m).ms_count = ar.model(m).ms_count + 1;
ar.model(m).data(d).ms_index = ar.model(m).ms_count;

for j=1:nints
    ar.model(m).data(d).ms_snip_index = j;
    if(j<nints)
        ar.model(m).data(end+1) = ar.model(m).data(d);
        ar.model(m).plot(jplot).dLink(end+1) = d+1;
    end
    
    if(j>1)
        ar.ms_count_snips = ar.ms_count_snips + 1;       
        qtodo = ismember(ar.model(m).data(d).p, ar.model(m).px0); %R2013a compatible
        ar.model(m).data(d).fp(qtodo) = strrep(ar.model(m).data(d).p(qtodo), 'init_', sprintf('init_MS%i_', ar.ms_count_snips));
    end
    
    if(j<nints)
        ar.model(m).data(d).tExp = ar.model(m).data(d).tExp(tExp>=tboarders(j) & tExp<tboarders(j+1));
        ar.model(m).data(d).yExp = ar.model(m).data(d).yExp(tExp>=tboarders(j) & tExp<tboarders(j+1),:);
        ar.model(m).data(d).yExpStd = ar.model(m).data(d).yExpStd(tExp>=tboarders(j) & tExp<tboarders(j+1),:);
    else
        ar.model(m).data(d).tExp = ar.model(m).data(d).tExp(tExp>=tboarders(j) & tExp<=tboarders(j+1));
        ar.model(m).data(d).yExp = ar.model(m).data(d).yExp(tExp>=tboarders(j) & tExp<=tboarders(j+1),:);
        ar.model(m).data(d).yExpStd = ar.model(m).data(d).yExpStd(tExp>=tboarders(j) & tExp<=tboarders(j+1),:);
    end
    
    ar.model(m).data(d).tLim = [tboarders(j) tboarders(j+1)];
    ar.model(m).data(d).tLimExp = ar.model(m).data(d).tLim;
    
    if(j<nints)
        d = d + 1;
    end
end


function ar = setValues(ar, m, d, header, nfactor, data, times)
ar.model(m).data(d).tExp = times;
ar.model(m).data(d).yExp = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStd = nan(length(times), length(ar.model(m).data(d).y));

for j=1:length(ar.model(m).data(d).y)
    q = ismember(header, ar.model(m).data(d).y{j}); %R2013a compatible
    
    if(sum(q)==1)
        ar.model(m).data(d).yExp(:,j) = data(:,q);
        fprintf('\t%20s -> %4i data-points assigned', ar.model(m).data(d).y{j}, sum(~isnan(data(:,q))));
        
        % normalize data
        if(ar.model(m).data(d).normalize(j))
            ar.model(m).data(d).yExp(:,j) = ar.model(m).data(d).yExp(:,j) / nfactor(q);
            fprintf(' normalized');
        end
        
        % log-fitting
        if(ar.model(m).data(d).logfitting(j))
            qdatapos = ar.model(m).data(d).yExp(:,j)>0;
            ar.model(m).data(d).yExp(qdatapos,j) = log10(ar.model(m).data(d).yExp(qdatapos,j));
            ar.model(m).data(d).yExp(~qdatapos,j) = nan;
            if(sum(~qdatapos)==0)
                fprintf(' for log-fitting');
            else
                fprintf(' for log-fitting (%i values <=0 removed)', sum(~qdatapos));
            end
        end
        
        % empirical stds
        qstd = ismember(header, [ar.model(m).data(d).y{j} '_std']); %R2013a compatible
        if(sum(qstd)==1)
            ar.model(m).data(d).yExpStd(:,j) = data(:,qstd);
            fprintf(' with stds');
            if(ar.model(m).data(d).normalize(j))
                ar.model(m).data(d).yExpStd(:,j) = ar.model(m).data(d).yExpStd(:,j) / nfactor(q);
                fprintf(' normalized');
            end
        elseif(sum(qstd)>1)
            error('multiple std colums for observable %s', ar.model(m).data(d).y{j})
        end
        
    elseif(sum(q)==0)
        fprintf('*\t%20s -> not assigned', ar.model(m).data(d).y{j});
    else
        error('multiple data colums for observable %s', ar.model(m).data(d).y{j})
    end
    
    fprintf('\n');
end

