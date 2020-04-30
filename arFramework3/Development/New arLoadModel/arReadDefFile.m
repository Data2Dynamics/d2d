
% D = arReadDefFile(name, [opts], [comment_string])
%
% This function combines the data structures created by user and read from file
%
%   name             char      filename of data definition file
%   comment_string   char      it is read from opts.
%                              ['//']
%   opts             strucr    to define data path
%                              ['Data/']


function D = arReadDefFile(name, opts, comment_string)

matVer = ver('MATLAB');

if ~exist('opts','var') || isempty(opts) ||isempty(opts.datapath_args)
    DataPath = 'Data/';
    opts={};
else
    DataPath = opts.datapath_args;
    if DataPath(end)~='/' && DataPath(end)~='\'
        DataPath = [DataPath,'/'];
    end
end

if ~exist('comment_string','var') || isempty(comment_string) || ~ischar(comment_string)
    arFprintf(3, '\nusing defaul comment string ''//'' \n');
    comment_string='//';
end


if ~exist('name','var') || isempty(name)
    error('Please provide a file name')
elseif ~ischar(name)
    error('File name should be string');
end

% load data def
if ~exist(DataPath,'dir')
    error('folder %s does not exist',DataPath)
end
if strcmp(strrep(name,' ',''),name)~=1
    disp(name)
    error('File names should not contain empty spaces. Please remove it.');
end


if(~exist([DataPath, name, '.def'],'file'))
    arFprintf(1, '\ncreating generic def file for %s%s ...\n', DataPath, name);
    D = arCreateDataStructure({},[],name,opts);
else
    
    D.name=name;
    D.path=[pwd '/' DataPath];
    D.condition=[];
    D.fp={};
    D.p={};
    D.py={};
    D.py_sep=struct;
    D.pystd = {};
    D.tExp = [];
    D.yExp = [];
    D.yExpRaw = [];
    D.yExpStd = [];
    D.yExpStdRaw = [];
    
    
    arFprintf(1, '\nLoading data definition from %s%s.def...', DataPath, name);
    
    % Disable this if you are having problems because of the preprocessor
    preprocessor = 1;
    if  ~preprocessor
        fid = fopen([DataPath, name, '.def'], 'r');
    else
        % Load into a struct
        fid.fn  = [DataPath name '.def'];
        fid.str = fileread([DataPath name '.def']);
        fid.pos = 1;
        arFprintf(3, 'Running preprocessor...' );
        fid = arPreProcessor(fid);
        arFprintf(3, ' [ OK ]\n' );
    end
    
    
    % DESCRIPTION ---------------------------------------------------------
    arFprintf( 3, 'Start parsing...' );
    [str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', comment_string);
    if ~strcmp(str{1},'DESCRIPTION')
        arParsingError(fid, 'parsing data %s for DESCRIPTION', name);
    end
    
    % check version
    if strcmp(str{1},'DESCRIPTION')
        % def_version = 1;
    elseif strcmp(str{1},'DESCRIPTION-V2')
        arParsingError(fid, 'DESCRIPTION-V2 not supported yet');
    else
        arParsingError(fid, 'invalid version identifier: %s', cell2mat(str{1}));
    end
    
    % read comments
    arFprintf( 3, '[ OK ]\nReading description' );
    [str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', comment_string);
    D.description = {};
    while ~strcmp(str{1},'PREDICTOR') && ~strcmp(str{1},'PREDICTOR-DOSERESPONSE')
        arFprintf(3, '.' );
        D.description(end+1,1) = str{1};
        [str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', comment_string);
    end
    
    
    % PREDICTOR -----------------------------------------------------------
    arFprintf( 3, '[ OK ]\nReading PREDICTOR ...\n' );
    if strcmp(str{1},'PREDICTOR-DOSERESPONSE')
        D.doseresponse = true;
        [str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', comment_string);
        D.response_parameter = cell2mat(str{1});
        arFprintf(2, 'dose-response to %s\n', D.response_parameter);
    else
        D.doseresponse = false;
        D.response_parameter = '';
        arFprintf(2, '\n');
    end
    [C, fid] = arTextScan(fid, '%s %s %q %q %n %n %n %n\n',1, 'CommentStyle', comment_string);
    
    D.t = cell2mat(C{1});
    D.tUnits(1) = C{2};
    D.tUnits(2) = C{3};
    D.tUnits(3) = C{4};
    D.tLim = [checkNum(C{5}, 0) checkNum(C{6}, 10)];
    D.tLimExp = [checkNum(C{7}, NaN), checkNum(C{8}, NaN)];
    
    
    % INPUTS --------------------------------------------------------------
    arFprintf( 3, '[OK]\nReading INPUTS ...\n' );
    [str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', comment_string);
    if ~strcmp(str{1},'INPUTS')
        arParsingError( fid, 'parsing data %s for INPUTS', name);
    end
    
    D.u = {};
    D.fu = {};
    D.uNames = {};
    D.pu={};
    [C, fid] = arTextScan(fid, '%s %q %q \n',1, 'CommentStyle', comment_string);
    while ~strcmp(C{1},'OBSERVABLES')
        
        arValidateInput(C, 'input', 'unique input name', 'input function', 'input label');
        D.u(end+1) = C{1};
        D.fu(end+1,1) = C{2};
        if ~isempty(cell2mat(C{3}))
            D.uNames(end+1) = C{3};
        else
            D.uNames{end+1} = '';
        end
        
        [C, fid] = arTextScan(fid, '%s %q %q \n',1, 'CommentStyle', comment_string);
    end
    
    
    
    % OBSERVABLES ---------------------------------------------------------
    arFprintf(3, '[ OK ]\nReading OBSERVALBES ...\n' );
    
    D.y = {};
    D.yNames = {};
    D.yUnits = {};
    D.normalize = [];
    D.logfitting = [];
    D.logplotting = [];
    D.fy = {};
    
    if strcmp(C{1},'OBSERVABLES')
        [C, fid] = arTextScan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', comment_string);
        while ~strcmp(C{1},'ERRORS')
            if strcmp(C{1}, 'CONDITIONS')
                arParsingError(fid, 'When OBSERVABLES section is specified; ERRORS section must also be specified.');
            end
            arValidateInput(C, 'observable', 'unique identifier', 'unit type (i.e. C)', 'unit (i.e. "units/cell")', 'plain text label for plots ("conc.")', 'indicator whether data should be scaled to 1 (0 or 1)', 'Indicator whether date should be treated in log-space (0 or 1)', 'Mathematical expression for observable' );
            D.y(end+1) = C{1};
            D.yUnits(end+1,1) = C{2};
            D.yUnits(end,2) = C{3};
            D.yUnits(end,3) = C{4};
            D.normalize(end+1) = C{5};
            D.logfitting(end+1) = C{6};
            D.logplotting(end+1) = C{6};
            D.fy(end+1,1) = C{7};
            if ~isempty(cell2mat(C{8}))
                D.yNames(end+1) = C{8};
            else
                D.yNames(end+1) = D.y(end);
            end
            
            if sum(ismember({D.t}, {D.y{end}}))>0 %R2013a compatible
                arParsingError(fid, '%s already defined in PREDICTOR', D.y{end});
            end
            if sum(ismember(D.u, D.y{end}))>0 %R2013a compatible
                arParsingError(fid, '%s already defined in INPUTS', D.y{end});
            end
            
            [C, fid] = arTextScan(fid, '%s %q %q %q %n %n %q %q\n',1, 'CommentStyle', comment_string);
        end
    end
    
    
    % ERRORS --------------------------------------------------------------
    arFprintf(3, '[ OK ]\nReading ERRORS ...\n');
    D.fystd={};
    if strcmp(C{1},'ERRORS')
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string);
        % put an error msg when observable is not empy but errors is empty
        while ~strcmp(C{1},'CONDITIONS') && ~strcmp(C{1},'DERIVED') && ~strcmp(C{1},'INVARIANTS') && ~strcmp(C{1},'SUBSTITUTIONS')
            arValidateInput(C, 'error', 'observable identifier', 'expression for the error model');
            qy = ismember(D.y, C{1}); %R2013a compatible
            if sum(qy)<1
                arParsingError(fid, 'Unknown observable in error specification %s.', cell2mat(C{1}));
            elseif sum(qy)>1
                arParsingError(fid, 'Multiple matches for %s', cell2mat(C{1}));
            end
            
            D.fystd(qy) = C{2};
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string);
        end
    end
    
    D.fystd=transpose(D.fystd);
    
    
    
    % DERIVED -------------------------------------------------------------
    if strcmp(C{1},'DERIVED')
        arParsingError(fid, ['There is no need for a section DERIVED in data definition file! ' ...
            'Please remove and see usage in: ' ...
            'https://github.com/Data2Dynamics/d2d/wiki/Setting%20up%20models']);
    end
    
    % INVARIANTS ----------------------------------------------------------
    if strcmp(C{1},'INVARIANTS')
        arParsingError( fid, ['Section INVARIANTS in data definition file is deprecated! ' ...
            'Please remove and see usage in: ' ...
            'https://github.com/Data2Dynamics/d2d/wiki/Setting%20up%20models']);
    end
    
    
    
    % SUBSTITUTIONS (beta) ------------------------------------------------
    
    if strcmp(C{1},'SUBSTITUTIONS')
        
        arFprintf(3, '[OK]\nReading SUBSTITUTIONS ...\n');
        if str2double(matVer.Version)>=8.4
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string, 'BufSize', 2^16);
        end
        
        % Substitutions
        D.subs.from = {};
        D.subs.to = {};
        
        % Fetch desired substitutions
        while ~isempty(C{1}) && ~strcmp(C{1},'CONDITIONS')
            arFprintf(3, '.' );
            D.subs.from(end+1) = C{1};
            D.subs.to(end+1) = C{2};
            
            if str2double(matVer.Version)>=8.4
                [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string);
            else
                [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string, 'BufSize', 2^16);
            end
            
        end
        
    end
    
    
    
    % CONDITIONS ----------------------------------------------------------
    D.pcond={};
    D.fpcond={};
    
    if strcmp(C{1},'CONDITIONS')
        
        arFprintf(3, '[OK]\nReading CONDITIONS ...\n');
        if str2double(matVer.Version)>=8.4
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string, 'BufSize', 2^16);
        end
        
        while ~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS') && ~strcmp(C{1},'RANDOM')
            arFprintf(3, '.' );
            arValidateInput(C, 'condition', 'model parameter', 'new expression');
            %qcondpara = ismember(arTemp.model(m).data(d).p, C{1}); %R2013a compatible
            D.pcond{end+1}  = cell2mat(C{1});
            D.fpcond{end+1} = ['(' cell2mat(C{2}) ')'];
            
            if str2double(matVer.Version)>=8.4
                [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string);
            else
                [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', comment_string, 'BufSize', 2^16);
            end
        end
        
    else
        arParsingError(fid, 'Missing field CONDITIONS');
    end
    
    D.pcond = transpose(D.pcond);
    D.fpcond= transpose(D.fpcond);
    
    
    % RANDOM --------------------------------------------------------------
    D.prand={};
    D.rand_type=[];
    
    if strcmp(C{1},'RANDOM')
        arFprintf(3, '[ OK ]\nReading RANDOM');
        [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', comment_string);
        
        while ~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS')
            arFprintf(3, '.' );
            D.prand{end+1} = cell2mat(C{1});
            if strcmp(C{2}, 'INDEPENDENT')
                D.rand_type(end+1) = 0;
            elseif strcmp(C{2}, 'NORMAL')
                D.rand_type(end+1) = 1;
            else
                warning('unknown random type %s', cell2mat(C{2}));  %#ok<WNTAG>
            end
            [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', comment_string);
        end
    end
    
    
    % PARAMETERS ----------------------------------------------------------
    D.par.pExternLabels = {};
    D.par.pExtern = [];
    D.par.qFitExtern = [];
    D.par.qLog10Extern = [];
    D.par.lbExtern = [];
    D.par.ubExtern = [];
    
    if strcmp(C{1},'PARAMETERS')
        arFprintf(3, '[ OK ]\nReading PARAMETERS...');
        [C, fid] = arTextScan(fid, '%s %f %n %n %n %n\n',1, 'CommentStyle', comment_string);
        while ~isempty(C{1})
            D.par.pExternLabels(end+1) = C{1};
            D.par.pExtern(end+1) = C{2};
            D.par.qFitExtern(end+1) = C{3};
            D.par.qLog10Extern(end+1) = C{4};
            D.par.lbExtern(end+1) = C{5};
            D.par.ubExtern(end+1) = C{6};
            [C, fid] = arTextScan(fid, '%s %f %n %n %n %n\n',1, 'CommentStyle', comment_string);
        end
        
    end
    
    
    % Close data def file
    if ~isstruct(fid)
        fclose(fid);
    end
    
end


D = orderfields(D);

end



function num = checkNum(num, defaultValue)
if ~isnumeric(num) || isempty(num) || isnan(num)
    num = defaultValue;
end
end