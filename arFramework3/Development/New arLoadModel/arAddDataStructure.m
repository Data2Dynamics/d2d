% arAddDataStructure(name, D, m, removeEmptyObs, opts, varargin)
% 
% This function updates and computes the data struct and appends it to ar.model(m). Compiling has
% to be performed after this step.


function arAddDataStructure(name, D, m, removeEmptyObs, opts, varargin)


global ar

matVer = ver('MATLAB');

if isempty(name) || endsWith(D.name,'AutoGen')
    name = D.name;
end

if ~exist('opts','var') || isempty(opts.datapath_args)
    DataPath = 'Data/';
else
    DataPath = opts.datapath_args;
    if DataPath(end)~='/' && DataPath(end)~='\'
        DataPath = [DataPath,'/'];
    end
end

if  ~exist('D','var') || isempty(D)
    error('please load data file')
end

if isempty(ar)
    error('please initialize by arInit')
end

if ~isfield(ar,'model')
    error('please load the model')
end


if ~exist('m','var') || isempty(m)
    m = length(ar.model);
end

if exist('m','var') && ischar(m)
    for jm=1:length(ar.model)
        if(strcmp(m, ar.model(jm).name))
            m = jm;
        end
    end
    if(ischar(m))
        error('Model %s was not found', m);
    end
end

% find/make an appropriate data field
if isfield(ar.model(m), 'data') 
    for i=length(ar.model(m).data)
        if strcmp(ar.model(m).data(i).name,D.name)
            d = i;
        else
            d = length(ar.model(m).data) + 1;
        end
    end
else
    %ar.model(m).data = {};
    d = 1;
end



if(ischar(removeEmptyObs))
    error(['arAddDataStruct(D, m ...) input argument removeEmptyObs is deprecated !!! ' ...
        'Please see new usage arAddDataStructure(name, DataDefInfo, m, removeEmptyObs, opts, varargin) and function help text.']);
end


if ( opts.resampledoseresponse )
    if ( ~isnumeric( opts.resamplingresolution_args ) || isempty( opts.resamplingresolution_args ) )
        opts.resamplingresolution = 25;
    else
        opts.resamplingresolution = opts.resamplingresolution_args(1);
    end
end

if( opts.dppershoot )
    if( opts.dppershoot_args>0 )
        if(~isfield(ar,'ms_count_snips'))
            ar.model(m).ms_count = 0;
            ar.ms_count_snips = 0;
            ar.ms_strength = 0;
            ar.ms_threshold = 1e-5;
            ar.ms_violation = [];
        end
        dpPerShoot = opts.dppershoot_args;
    end
else
    dpPerShoot = 0;
end



% PREDICTOR ---------------------------------------------------------------
D.tLimExp = [checkNum(D.tLimExp(1), ar.model(m).tLim(1)), checkNum(D.tLimExp(2), ar.model(m).tLim(2))];


% INPUTS ------------------------------------------------------------------
Dtemp.fu = ar.model.fu;
Dtemp.uNames = {};

for i=1:length(D.u)
    qu = ismember(Dtemp.u, D.u(i)); %R2013a compatible
    if(sum(qu)~=1)
        error('unknown input %s in file %s.def', D.u{i}, D.name );
    end    
    
    % Ignore this replacement?
    ignoreInput = 0;
    if (opts.ignoreinputs)
        if ismember(D.u(i), opts.ignoreinputs_args)
            ignoreInput = 1;
        end
    end
    if ( ~ignoreInput )
        % Input replacement description
        Dtemp.fu(qu) = D.fu(i);
        if(~isempty(D.uNames(i)))
            Dtemp.uNames(end+1) = D.uNames(i);
        else
            Dtemp.uNames{end+1} = '';
        end
    end    
end

D.fu = Dtemp.fu;
D.uNames = Dtemp.uNames;


% input parameters
varlist = cellfun(@symvar, D.fu, 'UniformOutput', false);
D.pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}); %R2013a compatible


% OBSERVABLES -------------------------------------------------------------
qy = ismember(ar.model(m).y,D.y);
idx=find(qy==0);

Dtemp.y=D.y;
Dtemp.fystd = D.fystd;

D.y = cat(2, D.y, ar.model(m).y(idx));
D.yNames = cat(2, D.yNames, ar.model(m).yNames(idx));
D.yUnits = cat(1, D.yUnits, ar.model(m).yUnits(idx,:));
D.normalize = cat(2, D.normalize, ar.model(m).normalize(idx));
D.logfitting = cat(2, D.logfitting, ar.model(m).logfitting(idx));
D.logplotting = cat(2, D.logplotting, ar.model(m).logplotting(idx));
D.fy = cat(1, D.fy, ar.model(m).fy(idx));

%error model
D.fystd = cat(1, D.fystd, ar.model(m).fystd(idx,:));


% we can do it using intersect(a,b)
if any(ismember(D.y, ar.model(m).x))
    error('%s already defined in STATES\n', D.y{find(ismember(D.y,ar.model(m).x)==1)});
end
if any(ismember(D.y, ar.model(m).u))
    error('%s already defined in INPUTS\n', D.y{find(ismember(D.y,ar.model(m).u)==1)});
end
if any(ismember(D.y, ar.model(m).z))
    error('%s already defined in DERIVED\n', D.y{find(ismember(D.y,ar.model(m).z)==1)});
end
if any(ismember(D.y, ar.model(m).p))
    error('%s already defined as parameter\n', D.y{find(ismember(D.y,ar.model(m).p)==1)});
end



% observation parameters
arFprintf( 3, '[ OK ]\nComputing observation parameters...' );
varlist = cellfun(@symvar, D.fy, 'UniformOutput', false);
D.py = setdiff(setdiff(vertcat(varlist{:}), union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
if(isempty(D.fy))
    error('No OBSERVABLE specified. Specify an OBSERVABLE in the model or data definition file. See "Defining the OBSERVABLES".');
end
for j=1:length(D.fy)
    varlist = symvar(D.fy{j});
    D.py_sep(j).pars = setdiff(setdiff(varlist, union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
    
    % exclude parameters form model definition
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).px);
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).pu);
end



% ERRORS ------------------------------------------------------------------

% check and error if observable in log and fystd = rel + abs error
for i=1:length(Dtemp.y)
    
    y_var_name = setdiff(symvar(D.fy{i}),D.py);
    reg_string = ['((?<=\W)|^)(',Dtemp.y{i},'|'];
    for jreg = 1:length(y_var_name)
        if(jreg<length(y_var_name))
            reg_string = [reg_string ,y_var_name{jreg},'|'];
        else
            reg_string = [reg_string ,y_var_name{jreg},')'];
        end
    end

    reg_string = [reg_string '((?=\W)|$)'];
    if(~isempty(regexp(Dtemp.fystd{i},reg_string,'ONCE')) && D.logfitting(i))
        warning(['You are trying to set up a relative error model within a log transformation. \n%s' ...
            'Comment out this error if you want to proceed anyway. To implement an absolute error in log, \n' ...
            'you can try the approach: \nyObs = sd_yObs + 1/2 * (a+sqrt((a)^2)), a = (offset - yObs-sd_yObs) \n, with hard set or fitted offset (on log-scale) \n'],C{2}{1})
        error('Revise error model')
    end
end


% error model D.fystd has been defined in OBSERVABLE section above!
if (length(D.fystd)<length(D.fy) || sum(cellfun(@isempty, D.fystd))>0)
    diffErr = D.y(cellfun(@isempty, D.fystd)>0);
    if ( length(D.fystd)<length(D.fy) )
        diffErr = union( D.y( length(D.fystd) + 1 : end ), diffErr );
    end
    error('Some observables do not have an error model defined. Observable(s) without error model: %s\n', sprintf( '%s ', diffErr{:} ) );
end


% Drop certain observables in OBSERVABLES & ERRORS 
if (opts.removeobservables)
    arFprintf( 3, '[ OK ]\nDropping specific observables...\n' );
    if ischar( opts.removeobservables_args )
        opts.removeobservables_args = {opts.removeobservables_args};
    end
    for a = 1 : length( opts.removeobservables_args )
        j = 1;
        while( j <= length( D.y ) )
            jind = ismember( D.y{j}, opts.removeobservables_args );
            if ( sum(jind) > 0 )
                warning( '>> Explicitly removing %s!\n', D.y{j} );
                D.y(j) = [];
                D.yUnits(j,:) = [];
                D.normalize(j) = [];
                D.logfitting(j) = [];
                D.logplotting(j) = [];
                D.fy(j) = [];
                D.yNames(j) = [];
                D.fystd(j) = [];
            else
                j = j + 1;
            end
        end
    end
end

% error parameters
arFprintf( 3, 'Compute error parameters...' );
varlist = cellfun(@symvar, D.fystd, 'UniformOutput', false);
D.pystd = setdiff(vertcat(varlist{:}), union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ... %R2013a compatible
    D.y), D.t));
for j=1:length(D.fystd)
    varlist = symvar(D.fystd{j});
	D.py_sep(j).pars = union(D.py_sep(j).pars, ... %R2013a compatible
        setdiff(varlist, union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), D.y))); %R2013a compatible
    
    % exclude parameters form model definition
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).px);
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).pu);
end


% collect parameters needed for OBSERVABLES
ptmp = union(ar.model(m).px, ar.model(m).pu);
D.p = union(ptmp, union(D.pu, D.py)); %R2013a compatible
D.pystd = setdiff(D.pystd, D.p); %Remove dynamic variables from error model parameters, e.g. scaling parameters if sd propto scale*x.
D.p = union(D.p, D.pystd); %R2013a compatible

% Union's behaviour is different when first arg is empty. In this case, a
% flip of the parameter vector is typically required.
if ( size( D.p, 1 ) ~= 1 )
    D.p = transpose(D.p); %instead of ar.model.data.p.'
end

% replace filename
D.p = strrep(D.p, '_filename', ['_' D.name]);
D.fy = strrep(D.fy, '_filename', ['_' D.name]);
D.py = strrep(D.py, '_filename', ['_' D.name]);
D.fystd = strrep(D.fystd, '_filename', ['_' D.name]);
D.pystd = strrep(D.pystd, '_filename', ['_' D.name]);
for j=1:length(D.py_sep)
    D.py_sep(j).pars = strrep(D.py_sep(j).pars, '_filename', ['_' D.name]);
end

clearvars Dtemp

% SUBSTITUTIONS (beta) ----------------------------------------------------
substitutions = 0;
if isfield(D,'subs') && ~isempty(D.subs.from)
    substitutions = 0;
    ismodelpar = ismember( D.subs.from, ar.model(m).p);
    
    if any(ismodelpar)
        error('Cannot substitute model parameter %s.\n', D.subs.from{ismodelpar} );
    end
    
    % Perform selfsubstitutions
    arFprintf( 3, '[ OK ]\nPerforming self substitutions...' );
    if ( ~isempty(D.subs.from) )
        substitutions = 1;
        D.subs.to = arSubsRepeated( D.subs.to, D.subs.from, D.subs.to, str2double(matVer.Version) );
    end
    arFprintf( 3, '[ OK ]\n' );
end



% CONDITIONS --------------------------------------------------------------
D.fp = transpose(D.p);
ptmp = ar.model(m).p;

qcondparamodel = ismember(D.p, strrep(ptmp, '_filename', ['_' D.name])); %R2013a compatible
qmodelparacond = ismember(strrep(ptmp, '_filename', ['_' D.name]), D.p); %R2013a compatible
D.fp(qcondparamodel) = strrep(ar.model(m).fp(qmodelparacond), '_filename', ['_' D.name]);

if ( substitutions == 1 )
    % Perform selfsubstitutions
    D.fpcond = arSubsRepeated( D.fpcond, D.subs.from, D.subs.to, str2double(matVer.Version) );
    
    % Store substitutions in ar structure
    for a = 1 : length( D.pcond )
        qcondpara = ismember(D.p, D.pcond(a)); %R2013a compatible
        if(sum(qcondpara)>0)
            D.fp{qcondpara} = ['(' D.fpcond{a} ')'];
        else
            warning('unknown parameter in conditions: %s (did you mean to place it under SUBSTITUTIONS?)', D.pcond{a}); %#ok<WNTAG>
        end
    end
    
else
    
    % old code path
    for a=1:length( D.pcond )
        arFprintf( 3, '.' );
        qcondpara = ismember(D.p, D.pcond(a)); %R2013a compatible
        if(sum(qcondpara)>0)
            D.fp{qcondpara} = ['(' cell2mat(D.fpcond(a)) ')'];
        else
            warning('unknown parameter in conditions %s', cell2mat(D.pcond{1}));
        end
    end    
    
end

% extra conditional parameters
varlist = cellfun(@symvar, D.fp, 'UniformOutput', false);
D.pcond = setdiff(vertcat(varlist{:}), D.p); %R2013a compatible
      
% collect parameters conditions
pcond = union(D.p, D.pcond); %R2013a compatible


% RANDOM --------------------------------------------------------------
qrand = ismember(ar.model(m).prand, D.prand);
idx=find(qrand==0);
D.prand = cat(2, D.prand, ar.model(m).prand(idx));
D.rand_type = cat(2, D.rand_type, ar.model(m).rand_type(idx));

if ( opts.expsplit )
    D.prand{end+1} = opts.expsplit_args;
    D.rand_type(end+1) = 0;
end

% PARAMETERS --------------------------------------------------------------
if(~isfield(ar, 'pExternLabels'))
    ar.pExternLabels = {};
    ar.pExtern = [];
    ar.qFitExtern = [];
    ar.qLog10Extern = [];
    ar.lbExtern = [];
    ar.ubExtern = [];
end

if isfield(D,'par')
    for i=1:length(D.par.pExternLabels)
        ar.pExternLabels(end+1) = D.par.pExternLabels(i);
        ar.pExtern(end+1) = D.par.pExtern(i);
        ar.qFitExtern(end+1) = D.par.qFitExtern(i);
        ar.qLog10Extern(end+1) = D.par.qLog10Extern(i);
        ar.lbExtern(end+1) = D.par.lbExtern(i);
        ar.ubExtern(end+1) = D.par.ubExtern(i);
    end
    D=rmfield(D,'par');
end


% plot setup --------------------------------------------------------------
arFprintf( 3, '[ OK ]\nPlot setup...' );
if(isfield(D, 'response_parameter') && ...
        ~isempty(D.response_parameter))
    if(sum(ismember(D.p ,D.response_parameter))==0 && ... %R2013a compatible
            sum(ismember(D.pcond ,D.response_parameter))==0) %R2013a compatible
        error('invalid response parameter %s', D.response_parameter);
    end
end
if(~isfield(ar.model(m), 'plot'))
    ar.model(m).plot(1).name = D.name;
else
    ar.model(m).plot(end+1).name = D.name;
end

ar.model(m).plot(end).doseresponse = D.doseresponse;
ar.model(m).plot(end).doseresponselog10xaxis = true;
ar.model(m).plot(end).dLink = d;
ar.model(m).plot(end).ny = length(D.y);
ar.model(m).plot(end).condition = {};
jplot = length(ar.model(m).plot);



% ar computation ----------------------------------------------------------
extension = D.extension;
data = D.xlsdata;
dataCell = D.xlsdataCell;
header = D.xlsheader;
times = D.xlstimes;
timevar = D.xlstimevar;

D = rmfield(D,'extension');
D = rmfield(D,'fpcond');
D = rmfield(D,'u');
D = rmfield(D,'xlsdata');
D = rmfield(D,'xlsdataCell');
D = rmfield(D,'xlsheader');
D = rmfield(D,'xlstimes');
D = rmfield(D,'xlstimevar');
%D.fprand='';

if ~isempty(m) && isfield(ar.model(m),'data') && ~isempty(ar.model(m).data)
    flds = fields(ar.model.data);    
    for i=1:length(flds)       
        if ~isfield(D,flds{i})         
            if iscell(ar.model.data.(flds{i}))
                D.(flds{i}) = cell(0);
            elseif isnumeric(ar.model.data.(flds{i}))
                D.(flds{i}) = [];
            elseif islogical(ar.model.data.(flds{i}))
                D.(flds{i}) = [];
            elseif ischar(ar.model.data.(flds{i}))
                D.(flds{i}) = '';
            elseif isstruct(ar.model.data.(flds{i}))
                D.(flds{i}) = struct;
            else
                error('This case is not yet implemented.');
            end
        end
    end     
end

ar.model(m).data(d) = D;


if ~strcmp(extension,'none')
    arFprintf(2, 'Computing data #%i ...\n', d);
    dataFound = true;
    
    % remove time points that we don't want
    if ( opts.removeconditions )
        arFprintf( 3, 'Removing undesired time points...' );
        selected = true(1, size(times,1));
        if ( opts.removeconditions )
            for a = 1 : 2 : length( opts.removeconditions_args )
                if ( strcmp( timevar, opts.removeconditions_args{a} ) )
                    % If the argument is a function handle, we evaluate them
                    % for each element
                    val = opts.removeconditions_args{a+1};
                    if ( isa(val, 'function_handle') )
                        for jv = 1 : length( times )
                            accepted(jv) = val(num2str(times(jv)));
                        end
                    else
                        error('Filter argument for removecondition is of the wrong type' );
                    end
                    selected = selected & ~accepted;
                end
            end
        end
        times    = times(selected);
        data     = data(selected,:);
        dataCell = dataCell(selected,:);
    end
    
       
    
    % random effects
    arFprintf( 3, 'Processing random effects...' );
    prand = ar.model(m).data(d).prand;
    if(opts.splitconditions)
        prand = union(prand, opts.splitconditions_args);
    end
    qrandis = ismember(header, prand); %R2013a compatible
    if(sum(qrandis) > 0)
        qobs = ismember(header, ar.model(m).data(d).y); %R2013a compatible
        
        randis_header = header(qrandis);
        qrandis_header_nosplit = ismember(randis_header, ar.model(m).data(d).prand);
        
        if ~isempty(dataCell)
            [randis, ~, jrandis] = uniqueRowsCA(dataCell(:,qrandis));
        else
            [randis, ~, jrandis] = unique(data(:,qrandis),'rows');
            randis = cellstr(num2str(randis));
        end
        
        for j=1:size(randis,1)
            qvals = jrandis == j;
            tmpdata = data(qvals,qobs);
            if(sum(~isnan(tmpdata(:)))>0 || ~removeEmptyObs)
                arFprintf(2, 'local random effect #%i:\n', j)
                
                if(j < size(randis,1))
                    ar.model(m).data(d+1) = ar.model(m).data(d);
                    ar.model(m).plot(jplot+1) = ar.model(m).plot(jplot);
                end
                
                pcondmod = pcond;
                for jj=1:size(randis,2)
                    if(qrandis_header_nosplit(jj))
                        arFprintf(2, '\t%20s = %s\n', randis_header{jj}, randis{j,jj})
                        
                        ar.model(m).plot(jplot).name = [ar.model(m).plot(jplot).name '_' ...
                            randis_header{jj} randis{j,jj}];
                         
                        ar.model(m).data(d).name = [ar.model(m).data(d).name '_' ...
                            randis_header{jj} randis{j,jj}];

                        ar.model(m).data(d).fprand = randis{j,jj};
                        
                        ar.model(m).data(d).fy = strrep(ar.model(m).data(d).fy, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).py = strrep(ar.model(m).data(d).py, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        
                        ar.model(m).data(d).fystd = strrep(ar.model(m).data(d).fystd, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).pystd = strrep(ar.model(m).data(d).pystd, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
%                         
%                         ar.model(m).data(d).p = strrep(ar.model(m).data(d).p, ...
%                             randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).fp = strrep(ar.model(m).data(d).fp, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        ar.model(m).data(d).pcond = strrep(ar.model(m).data(d).pcond, ...
                            randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        
                        for jjj=1:length(ar.model(m).data(d).py_sep)
                            ar.model(m).data(d).py_sep(jjj).pars = strrep(ar.model(m).data(d).py_sep(jjj).pars, ...
                                randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                        end
                        
                        pcondmod = strrep(pcondmod, randis_header{jj}, [randis_header{jj} randis{j,jj}]);
                    else
                        arFprintf(2, '\t%20s (split only)\n', randis_header{jj})
                        
                        ar.model(m).plot(jplot).name = [ar.model(m).plot(jplot).name '_' ...
                            randis_header{jj} randis{j,jj}];
                    end
                end
                
                arFprintf( 4, 'Setting conditions ...\n' );
                if ~isempty(dataCell)
                    [ar,d,fail] = setConditions(ar, m, d, jplot, header, ...
                        times(qvals), data(qvals,:), ...
                        dataCell(qvals,:), pcondmod, removeEmptyObs, dpPerShoot, opts);
                else
                    [ar,d,fail] = setConditions(ar, m, d, jplot, header, ...
                        times(qvals), data(qvals,:), ...
                        dataCell, pcondmod, removeEmptyObs, dpPerShoot, opts);
                end
                arFprintf( 4, 'Condition set ... [ OK ]\n' );
                
                % Only increment if some data was actually set.
                if (~fail)
                    if(j < size(randis,1))
                        d = d + 1;
                        jplot = jplot + 1;
                        ar.model(m).plot(jplot).dLink = d;
                    end
                    
                    % Check whether the user specified any variables with reserved words.
                    checkReserved(m, d);
                else
                    % Remove file which failed to provide any data
                    fprintf(2, 'local random effect #%i: no matching data (%d), removed\n', j ,d);
                    ar.model.data(d) = [];
                end
                
            else
                arFprintf(2, 'local random effect #%i: no matching data, skipped\n', j);
            end
        end
    else
        arFprintf( 4, 'Setting conditions ...\n' );
        ar = setConditions(ar, m, d, jplot, header, times, ...
            data, dataCell, pcond, removeEmptyObs, dpPerShoot, opts);
        arFprintf( 4, 'Condition set ... [ OK ]\n' );
        
        % Check whether the user specified any variables with reserved words.
        checkReserved(m, d);
    end
    
else
    dataFound = false;
    warning('Cannot find data file corresponding to %s', name);
    ar.model(m).data(d).condition = [];
end
    


% remember the function call
% ar.setup.commands{end+1} = mfilename; % this file name, this is done in
% new arLoadModel
ar.setup.arguments{end+1} = {ar.model(m).data(d).name,m,extension, removeEmptyObs, varargin{:}};
if dataFound
    ar.setup.datafiles{end+1} = {[DataPath name '.def'],[DataPath, name,'.',extension]};
else
    ar.setup.datafiles{end+1} = {[DataPath,name,'.def'],''};
end
ar.setup.modelfiles{end+1} = '';


% sort fields
ar = orderfields(ar);
ar.model = orderfields(ar.model);
ar.model(m).data = orderfields(ar.model(m).data);
ar.model(m).plot = orderfields(ar.model(m).plot);

end


function checkReserved(m, d)
    global ar;

    % Check whether the user specified any variables with reserved words.
    for a = 1 : length( ar.model(m).data(d).fu )
        arCheckReservedWords( symvar(ar.model(m).data(d).fu{a}), sprintf( 'input function of %s', ar.model(m).data(d).name ), ar.model(m).u{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fy )
        arCheckReservedWords( symvar(ar.model(m).data(d).fy{a}), sprintf( 'observation function of %s', ar.model(m).data(d).name ), ar.model(m).data(d).y{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fystd )
        arCheckReservedWords( symvar(ar.model(m).data(d).fystd{a}), sprintf( 'observation standard deviation function of %s', ar.model(m).data(d).name ), ar.model(m).data(d).y{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fp )
        arCheckReservedWords( symvar(ar.model(m).data(d).fp{a}), sprintf( 'condition parameter transformations of %s', ar.model(m).data(d).name ), ar.model(m).data(d).p{a} );
    end   
    arCheckReservedWords( ar.model(m).data(d).p, 'parameters' );
    arCheckReservedWords( ar.model(m).data(d).y, 'observable names' );
end

function num = checkNum( num, defaultValue )
    if ( ~isnumeric( num ) || isempty( num ) || isnan( num ) )
        num = defaultValue;
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

arFprintf(2, 'using %i shooting intervals\n', nints);
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
end

function C = mymat2cell(D)
C = cell(size(D));
for j=1:size(D,1)
    for jj=1:size(D,2)
        C{j,jj} = num2str(D(j,jj));
    end
end
end

function [ar,d, fail] = setConditions(ar, m, d, jplot, header, times, data, dataCell, pcond, removeEmptyObs, dpPerShoot, opts)

% normalization of columns
fail = 0;
nfactor = max(data, [], 1);

qobs = ismember(header, ar.model(m).data(d).y) & sum(~isnan(data),1)>0; %R2013a compatible
qhasdata = ismember(ar.model(m).data(d).y, header(qobs)); %R2013a compatible

% conditions
if (~opts.removeconditions)
    qcond = ismember(header, pcond); %R2013a compatible
else
    % Add the condi's we force filtering over (override)
    qcond = ismember(header, pcond) | ismember(header, opts.removeconditions_args(1:2:end)); %R2013a compatible
end

% Refine dose responses if requested
if ( opts.resampledoseresponse )
    resolution = opts.resamplingresolution;
    if ( ar.model(m).data(d).doseresponse == true )
        responsePar = ismember( header, ar.model(m).data(d).response_parameter );
        if ( sum( responsePar ) )
            fprintf( '  => Refining dose response for %s\n', ar.model(m).data(d).response_parameter );
            if ( sum( responsePar ) > 1 )
                error('Response parameter ambiguous during dose response refinement' );
            end

            % Which columns define the conditions
            conds = qcond & ~responsePar;

            % Grab unique conditions (note the inclusion of time)
            [uniqueCondi, ~, ib] = unique( [times data( :, qcond & ~responsePar ) ], 'rows' );
            nConditions = size(uniqueCondi, 1);

            % For each unique condition determine the maximum and minimum value
            % of the response parameter
            extraData = NaN(nConditions*resolution, size(data,2));
            extraTimes = NaN(nConditions*resolution, 1);

            for jui = 1 : size( uniqueCondi, 1 )
                dataChunk = data( ib == jui, responsePar );
                mi = min( dataChunk );
                ma = max( dataChunk );
                extraPoints = mi : (ma-mi)/(resolution-1) : ma;

                if ( opts.refinelog )     
                    mi = log10( max( [mi, 1e-8] ) );
                    ma = log10( max( [ma, 1e-8] ) );
                    extraPoints = 10.^(mi : (ma-mi)/(resolution-1) : ma);
                end
                
                % Fix to make sure the length is correct for the next
                % assignment. This is ok, since the unique will remove
                % these points again at a later stage.
                if ( mi == ma )
                    extraPoints = repmat( mi, 1, resolution );
                end

                % Fill extra data with current condition
                condition = uniqueCondi(jui,:);
                extraData((jui-1)*resolution+1:jui*resolution, conds) = repmat(condition(2:end), resolution, 1);
                extraTimes((jui-1)*resolution+1:jui*resolution) = repmat(condition(1), resolution, 1);

                % Fill the dependent variable with the new dependent variable values
                extraData((jui-1)*resolution+1:jui*resolution, responsePar) = extraPoints;
            end
        end
        data = [ data ; extraData ];
        dataCell = [ dataCell; cellfun(@num2str,num2cell(extraData), 'UniformOutput', false) ];
        times = [ times ; extraTimes ];
    end
end

if(sum(qcond) > 0)
    condi_header = header(qcond);
    if ~isempty(dataCell)
        [condis, ind, jcondis] = uniqueRowsCA(dataCell(:,qcond));
    else
        [condis, ind, jcondis] = unique(data(:,qcond),'rows');
        condis = mymat2cell(condis);
    end

    if (opts.removeconditions || opts.removeemptyconds)
        selected = true(1, size(condis,1));
        if ( opts.removeconditions )
            for a = 1 : 2 : length( opts.removeconditions_args )
                cc = ismember( condi_header, opts.removeconditions_args{a} );
                if ( sum( cc ) > 0 )
                    values = condis(:,cc);

                    % If the argument is a function handle, we evaluate them
                    % for each element
                    val = opts.removeconditions_args{a+1};
                    if ( isa(val, 'function_handle') )
                        for jv = 1 : length( values )
                            accepted(jv) = val(values{jv});
                        end
                    else
                        if (isnumeric(val))
                            val = num2str(val);
                        end
                        if ~ischar(val)
                            error('Filter argument for removecondition is of the wrong type' );
                        end
                        accepted = ismember(values, val).';
                    end
                    selected = selected & ~accepted;
                end
            end
        end
        if(opts.removeemptyconds)
            % Find out for which conditions we actually have data
            hasD = max(~isnan(data(ind,qobs)), [], 2);
            selected(hasD==0) = false;
        end
        condis = condis(selected,:);
        
        % Recompute jcondi's (list which points which data row corresponds
        % to which condition.
        mapTo   = cumsum(selected);
        mapTo(~selected) = -1;
        jcondis = mapTo(jcondis);
    end
    
    % exit if no data left
    if(size(condis,1)==0)
        fail = 1;
        return
    end
       
    active_condi = false(size(condis(1,:)));
    tmpcondi = condis(1,:);
    for j1=2:size(condis,1)
        for j2=1:size(condis,2)
            active_condi(j2) = active_condi(j2) | (~strcmp(tmpcondi{j2}, condis{j1,j2}));
        end
    end
        
    for j=1:size(condis,1)
        
        arFprintf(2, 'local condition #%i:\n', j)
        
        if(j < size(condis,1))
            if(length(ar.model(m).data) > d)
                ar.model(m).data(d+2) = ar.model(m).data(d+1);
            end
            ar.model(m).data(d+1) = ar.model(m).data(d);
        end
        
        % remove obs without data
        if(removeEmptyObs)
            for jj=find(~qhasdata)
                arFprintf(2, '\t%20s no data, removed\n', ar.model(m).data(d).y{jj});
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
                arFprintf(2, '\t%20s = %s\n', condi_header{jj}, condis{j,jj})
                
                qcondjj = ismember(ar.model(m).data(d).p, condi_header{jj}); %R2013a compatible
                if(sum(qcondjj)>0)
                    ar.model(m).data(d).fp{qcondjj} =  ['(' condis{j,jj} ')'];
                end
                qcondjj = ~strcmp(ar.model(m).data(d).p, ar.model(m).data(d).fp');
                if(~isnan(str2double(condis{j,jj})))
%                     ar.model(m).data(d).fp(qcondjj) = strrep(ar.model(m).data(d).fp(qcondjj), ...
%                         condi_header{jj}, condis{j,jj});

                    ar.model(m).data(d).fp(qcondjj) = regexprep(ar.model(m).data(d).fp(qcondjj),...
                        sprintf('\\<%s\\>', condi_header{jj}),condis{j,jj});
                    
%                     tmpfp = subs(sym(ar.model(m).data(d).fp(qcondjj)), ...
%                         sym(condi_header{jj}), sym(condis{j,jj}));
%                     jps = find(qcondjj);
%                     for jp = 1:length(jps)
%                         ar.model(m).data(d).fp{jps(jp)} = char(tmpfp(jp));
%                     end
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
            arFprintf(2, '\t%20s no data, removed\n', ar.model(m).data(d).y{jj});
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
    ar.model(m).data(d).tLim(2) = round(max(times)*1.1);
    
    if(dpPerShoot~=0)
        [ar,d] = doMS(ar,m,d,jplot,dpPerShoot);
    end
end

end

function ar = setValues(ar, m, d, header, nfactor, data, times)
ar.model(m).data(d).tExp = times;
ar.model(m).data(d).yExp = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStd = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpRaw = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStdRaw = nan(length(times), length(ar.model(m).data(d).y));

for j=1:length(ar.model(m).data(d).y)
    q = ismember(header, ar.model(m).data(d).y{j}); %R2013a compatible
    
    if(sum(q)==1)
        ar.model(m).data(d).yExp(:,j) = data(:,q);
        ar.model(m).data(d).yExpRaw(:,j) = data(:,q);
        arFprintf(2, '\t%20s -> %4i data-points assigned', ar.model(m).data(d).y{j}, sum(~isnan(data(:,q))));
        
        % normalize data
        if(ar.model(m).data(d).normalize(j))
            ar.model(m).data(d).yExp(:,j) = ar.model(m).data(d).yExp(:,j) / nfactor(q);
            arFprintf(2, ' normalized');
        end
        
        % log-fitting
        if(ar.model(m).data(d).logfitting(j))
            qdatapos = ar.model(m).data(d).yExp(:,j)>0;
            nancount = length(isnan(ar.model(m).data(d).yExp(:,j)));
            ar.model(m).data(d).yExp(qdatapos,j) = log10(ar.model(m).data(d).yExp(qdatapos,j));
            ar.model(m).data(d).yExp(~qdatapos,j) = nan;
            if(sum(~qdatapos)==0)
                arFprintf(2, ' for log-fitting');
            else
                arFprintf(2, ' for log-fitting (%i values <=0 removed, %i NaN values removed)', sum(~qdatapos)-nancount, nancount);
            end
        end
        
        % empirical stds
        qstd = ismember(header, [ar.model(m).data(d).y{j} '_std']); %R2013a compatible
        if(sum(qstd)==1)
            ar.model(m).data(d).yExpStdRaw(:,j) = data(:,qstd);
            ar.model(m).data(d).yExpStd(:,j) = data(:,qstd);
            arFprintf(2, ' with stds');
            if(ar.model(m).data(d).normalize(j))
                ar.model(m).data(d).yExpStd(:,j) = ar.model(m).data(d).yExpStd(:,j) / nfactor(q);
                arFprintf(2, ' normalized');
            end
        elseif(sum(qstd)>1)
            error('multiple std colums for observable %s', ar.model(m).data(d).y{j})
        end
        
    elseif(sum(q)==0)
        arFprintf(2, '*\t%20s -> not assigned', ar.model(m).data(d).y{j});
    else
        error('multiple data colums for observable %s', ar.model(m).data(d).y{j})
    end
    
    arFprintf(1, '\n');
end
end