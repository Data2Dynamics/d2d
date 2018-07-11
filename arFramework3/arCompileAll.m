% Compile c- and mex-files for models, conditions and data sets
%
% arCompileAll(forcedCompile, debug_mode, ext_source_dir)
%   forcedCompile:                                      [false]
%   debug_mode:         exclude precompiled objects     [false]
%   source_dir:         source directory                []
%   simplifyEquations:  simplify equations              [0]
%                        0 - Do not simplify
%                        1 - Only simplify fz
%                        2 - Simplify fz and fv (slow!)
%
% See https://github.com/Data2Dynamics/d2d/wiki/First%20steps 
% for description of work flow. 
%
% Copyright Andreas Raue 2013 (andreas.raue@fdm.uni-freiburg.de)

function arCompileAll(forcedCompile, debug_mode, source_dir, simplifyEquations)

% turn off warning message for custom c functions

warnreset = warning;
warning('off', 'symbolic:mupadmex:MuPADTextWarning');
warning('off', 'symbolic:generate:FunctionNotVerifiedToBeValid')
warning('off', 'symbolic:sym:sym:DeprecateExpressions')

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

[~, c_version_code] = arGetVersion;
if(~strcmp(c_version_code, ar.info.c_version_code))
    error('Workspace c-code version (%s) does not match framework c-code version (%s), please rerun setup script or downgrade code version.', ...
        c_version_code, ar.info.c_version_code);
end
if(~exist('forcedCompile','var'))
    forcedCompile = false;
end
if(~exist('debug_mode','var'))
    debug_mode = false;
end
if(~exist('source_dir','var'))
    source_dir = cd;
else
    if ( isempty(source_dir) )
        source_dir = cd;
    end
end
if(~exist('simplifyEquations', 'var'))
    simplifyEquations = 0;
end
if (isfield(ar.config, 'legacy_steps'))
    legacy_steps = ar.config.legacy_steps;
else
    legacy_steps = 0;
end

% remember the function call
ar.setup.commands{end+1} = mfilename; % this file name
ar.setup.arguments{end+1} = {forcedCompile, debug_mode, source_dir, simplifyEquations}; % 
ar.setup.datafiles{end+1} = {};
ar.setup.modelfiles{end+1} = '';


if(isfield(ar,'checkstr'))
    prepareBecauseOfRepeatedCompilation;
end

% Special function definitions
if ( ~legacy_steps )
    % Argument formats for functional replacements.
    % First indicates the function name, second indicates the function,
    % third is the argument list mapping. Note that only valid symbolic
    % expressions can be used as functional replacements!
    % The smooth steps are based on fermi functions. Note that the location
    % parameters have to be in increasing order for this to work.
    %
    % Note that they have to be ordered in decreasing identifier length
    % (first column); since detection is based on string searches.
    ar.config.specialFunc = { ...
            {'smoothstep1',     '%s + (%s-%s) / (exp((%s-%s) / %s) + 1)', [4, 2, 4, 1, 3, 5], 'smoothstep1(t, level1, switch_time, level2, smoothness)' }, ...
            {'smoothstep2',     '%s + (%s-%s) / (exp((%s-%s) / %s) + 1) + (%s-%s) / (exp((%s-%s) / %s) + 1)', [6, 2, 4, 1, 3, 7, 4, 6, 1, 5, 7], 'smoothstep2( t, level1, switch_time1, level2, switch_time2, level3, smoothness )' }, ...        
            {'step1',           '%s + (%s-%s) * heaviside(%s-%s)', [2, 4, 2, 1, 3], 'step1(t, level1, switch_time, level2)'}, ...
            {'step2',           '%s + (%s-%s) * heaviside(%s-%s) + (%s-%s)*heaviside(%s-%s)', [2, 4, 2, 1, 3, 6, 4, 1, 5], 'step2(t, level1, switch_time1, level2, switch_time2, level3)' }, ...
            {'bolus',           '%s * (1 / sqrt( 2 * pi * %s^2 ) ) * exp(-(%s - %s)^2 / (2*%s^2))', [2, 4, 1, 3, 4], 'bolus(t, amount, time_point, duration)' }, ...
            {'hill_ka',         '%s^%s / (%s^%s + %s^%s)', [1, 3, 2, 3, 1, 3], 'hill_ka(conc, ka, n )' }, ...
            {'hill_kd',         '%s^%s / (%s + %s^%s)', [1, 3, 2, 1, 3], 'hill_kd(conc, kd, n )' }, ...
            {'isnonzero',       '(2*heaviside(%s))-1', 1, 'isnonzero(level)'}, ...
            {'max2',            '0.5*(%s+%s+abs(%s-%s))', [1,2,1,2], 'max2(a, b)'}, ...
            {'smooth1',         '(heaviside((%s - %s)/(%s - %s)) - (heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s))^2*((2*heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s) - 2*heaviside((%s - %s)/(%s - %s)) + 3)', [3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3], 'smooth1(t, start, end)' }, ...
            {'smooth2',         '-(heaviside((%s - %s)/(%s - %s)) - (heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s))^3*((heaviside((%s - %s)/(%s - %s)) - (heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s))*((6*heaviside((%s - %s)/(%s - %s))*(%s - %s)*(heaviside((%s - %s)/(%s - %s)) - 1))/(%s - %s) - 6*heaviside((%s - %s)/(%s - %s)) + 15) - 10)', [3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 1, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3], 'smooth2(t, start, end)' }, ...
        };
    
    % Add brackets for replacement safety
    for a = 1 : size( ar.config.specialFunc, 2 )
        ar.config.specialFunc{a}{2} = strrep(ar.config.specialFunc{a}{2}, '%s', '(%s)');
    end
else
    ar.config.specialFunc = [];
end

% folders
if(~exist([source_dir '/Compiled'], 'dir'))
	mkdir([source_dir '/Compiled'])
end
if(~exist([source_dir '/Compiled/' ar.info.c_version_code], 'dir'))
	mkdir([source_dir '/Compiled/' ar.info.c_version_code])
end

% Compiled folder hook for cluster usage
if(~exist([source_dir '/Compiled/arClusterCompiledHook.m'],'file'))
    fid = fopen([source_dir '/Compiled/arClusterCompiledHook.m'], 'W');
    fprintf(fid, 'function arClusterCompiledHook\n');
    fclose(fid);
end

% enable timedebug mode, use with debug_mode = true!
timedebug = false;

usePool = exist('gcp','file')>0 && ~isempty(gcp('nocreate'));

% main loop
checksum_global = addToCheckSum(ar.info.c_version_code);
c_version_code = ar.info.c_version_code;
for m=1:length(ar.model)
    arFprintf(2, '\n');
    
    % Have we reduced the model?
    if isfield( ar.model(m), 'reducedForm' )
        checksum_global = addToCheckSum(ar.model(m).reducedForm, checksum_global);
        
        % Store whether we express the model in totals or not
        for jt = 1 : numel( ar.model(m).removedStates )
            checksum_global = addToCheckSum(ar.model(m).removedStates(jt).expressInTotal, checksum_global);
        end
    end
    
    matVer = ver('MATLAB');
    matlab_version = str2double(matVer.Version);    
    
    % calc model
    arCalcModel(m, matlab_version, ar.config.networkgraph);
    
    % extract conditions
    ar.model(m).condition = [];
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            
            % conditions checksum
            qdynparas = ismember(ar.model(m).data(d).p, ar.model(m).px) | ... %R2013a compatible
                ismember(ar.model(m).data(d).p, ar.model(m).data(d).pu); %R2013a compatible
            
            checksum_cond = addToCheckSum(ar.model(m).data(d).fu);
            checksum_cond = addToCheckSum(ar.model(m).px, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).fv, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).N, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).cLink, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).z, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).fz, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).data(d).fp(qdynparas), checksum_cond);
            if isfield( ar.model(m), 'reducedForm' )
                % Store whether we express the model in totals or not
                for jt = 1 : numel( ar.model(m).removedStates )
                    checksum_cond = addToCheckSum(ar.model(m).removedStates(jt).expressInTotal, checksum_cond);
                end
            end
            
            checkstr_cond = getCheckStr(checksum_cond);
            
            % data checksum
            checksum_data = addToCheckSum(ar.model(m).data(d).fu);
            checksum_data = addToCheckSum(ar.model(m).data(d).p, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fy, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fystd, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fp, checksum_data);
            checkstr_data = getCheckStr(checksum_data);
            
            ar.model(m).data(d).checkstr = checkstr_data;
            ar.model(m).data(d).fkt = [ar.model(m).data(d).name '_' checkstr_data];
            
            cindex = -1;
            for c=1:length(ar.model(m).condition)
                if(strcmp(checkstr_cond, ar.model(m).condition(c).checkstr))
                    cindex = c;
                end
            end
            
            % global checksum
            if(isempty(checksum_global))
                checksum_global = addToCheckSum(ar.model(m).data(d).fkt);
            else
                checksum_global = addToCheckSum(ar.model(m).data(d).fkt, checksum_global);
            end
            
            if(cindex == -1) % append new condition
                cindex = length(ar.model(m).condition) + 1;
                
                ar.model(m).condition(cindex).status = 0;
                
                ar.model(m).condition(cindex).fu = ar.model(m).data(d).fu;
                ar.model(m).condition(cindex).fp = ar.model(m).data(d).fp(qdynparas);
                ar.model(m).condition(cindex).p = ar.model(m).data(d).p(qdynparas);
                
                ar.model(m).condition(cindex).checkstr = checkstr_cond;
                ar.model(m).condition(cindex).fkt = [ar.model(m).name '_' checkstr_cond];
                
                ar.model(m).condition(cindex).dLink = d;
                
                % global checksum
                checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt, checksum_global);
                
                % link data to condition
                ar.model(m).data(d).cLink = length(ar.model(m).condition);
                
                % for multiple shooting
                if(isfield(ar.model(m).data(d), 'ms_index') && ~isempty(ar.model(m).data(d).ms_index))
                    ar.model(m).condition(cindex).ms_index = ...
                        ar.model(m).data(d).ms_index;
                    ar.model(m).condition(cindex).ms_snip_index = ...
                        ar.model(m).data(d).ms_snip_index;
                    ar.model(m).condition(cindex).ms_snip_start = ar.model(m).data(d).tLim(1);
                end
            else
                % link data to condition
                ar.model(m).condition(cindex).dLink(end+1) = d;
                ar.model(m).data(d).cLink = cindex;
                
                % for multiple shooting
                if(isfield(ar.model(m).data(d), 'ms_index') && ~isempty(ar.model(m).data(d).ms_index))
                    ar.model(m).condition(cindex).ms_index(end+1) = ...
                        ar.model(m).data(d).ms_index;
                    ar.model(m).condition(cindex).ms_snip_index(end+1) = ...
                        ar.model(m).data(d).ms_snip_index;
                    ar.model(m).condition(cindex).ms_snip_start(end+1) = ar.model(m).data(d).tLim(1);
                end
            end
        end
        
        % skip calc conditions
        doskip = nan(1,length(ar.model(m).condition));
        for c=1:length(ar.model(m).condition)
            doskip(c) = ~forcedCompile && exist([source_dir '/Compiled/' ar.info.c_version_code '/' ar.model(m).condition(c).fkt '.c'],'file');
        end
        
        % calc conditions
        config = ar.config;
        model.name = ar.model(m).name;
        model.fv = ar.model(m).fv;
        model.fz = ar.model(m).fz;
        model.px0 = ar.model(m).px0;
        model.sym = ar.model(m).sym;
        model.t = ar.model(m).t;
        model.x = ar.model(m).x;
        model.u = ar.model(m).u;
        model.z = ar.model(m).z;
        model.us = ar.model(m).us;
        model.xs = ar.model(m).xs;
        model.vs = ar.model(m).vs;
        model.zs = ar.model(m).zs;
        model.N = ar.model(m).N;
        
        model.dvdx = ar.model(m).sym.dvdx;
        model.dvdu = ar.model(m).sym.dvdu;
        
        if isfield( ar.model(m), 'removedStates' )
            arFprintf( 2, 'Dealing with reduced model ...\n' );
            model.removedStates = ar.model(m).removedStates;
        end
        condition = ar.model(m).condition;
        newp = cell(1,length(ar.model(m).condition));
        newfp = cell(1,length(ar.model(m).condition));
        newpold = cell(1,length(ar.model(m).condition));
        newpx0 = cell(1,length(ar.model(m).condition));
        splines = cell(1,length(ar.model(m).condition));
        constVars = cell(1, length(ar.model(m).condition));
        
        if(usePool)
            csyms = cell(size(ar.model(m).condition));
            parfor c=1:length(ar.model(m).condition)
                warning('off','symbolic:mupadmex:MuPADTextWarning');
                warning('off', 'symbolic:generate:FunctionNotVerifiedToBeValid')
                warning('off', 'symbolic:sym:sym:DeprecateExpressions')
                
                condition_sym = arCalcCondition(config, model, condition(c), m, c, doskip(c), matlab_version, simplifyEquations);
                csyms{c} = condition_sym.sym;
                newp{c} = condition_sym.p;
                newfp{c} = condition_sym.fp;
                newpold{c} = condition_sym.pold;
                newpx0{c} = condition_sym.px0;
                splines{c} = condition_sym.splines;
                constVars{c} = condition_sym.constVars;
                
                if(~doskip(c))
                    % header
                    fid_odeH = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], 'W');
                    arWriteHFilesCondition(fid_odeH, config, condition_sym);
                    fclose(fid_odeH);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.h'],'f');
                    % body
                    fid_ode = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], 'W');
                    arWriteCFilesCondition(fid_ode, matlab_version, config, model, condition_sym, m, c, timedebug);
                    fclose(fid_ode);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.c'],'f');
                end
            end
            for c=1:length(condition)
                ar.model(m).condition(c).sym = csyms{c};
            end
        else
            for c=1:length(ar.model(m).condition)
                condition_sym = arCalcCondition(config, model, condition(c), m, c, doskip(c), matlab_version, simplifyEquations);
                ar.model(m).condition(c).sym = condition_sym.sym;
                newp{c} = condition_sym.p;
                newfp{c} = condition_sym.fp;
                newpold{c} = condition_sym.pold;
                newpx0{c} = condition_sym.px0;
                splines{c} = condition_sym.splines;
                constVars{c} = condition_sym.constVars;
                
                if(~doskip(c))
                    % header
                    fid_odeH = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], 'W');
                    arWriteHFilesCondition(fid_odeH, config, condition_sym);
                    fclose(fid_odeH);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.h'],'f');
                    % body
                    fid_ode = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], 'W');
                    arWriteCFilesCondition(fid_ode, matlab_version, config, model, condition_sym, m, c, timedebug);
                    fclose(fid_ode);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.c'],'f');
                end
            end
        end
        
        % assign conditions
        for c=1:length(ar.model(m).condition)
            ar.model(m).condition(c).p = newp{c};
            ar.model(m).condition(c).fp = newfp{c};
            ar.model(m).condition(c).pold = newpold{c};
            ar.model(m).condition(c).px0 = newpx0{c};
            ar.model(m).condition(c).splines = splines{c};
            ar.model(m).condition(c).constVars = constVars{c};
        end
        
        % skip calc data
        doskip = nan(1,length(ar.model(m).data));
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).p_condition = ar.model(m).condition(ar.model(m).data(d).cLink).p;
            doskip(d) = ~forcedCompile && exist([source_dir '/Compiled/' ar.info.c_version_code '/' ar.model(m).data(d).fkt '.c'],'file');
        end
        
        % calc data
        data = ar.model(m).data;
        newp = cell(1,length(ar.model(m).data));
        newfp = cell(1,length(ar.model(m).data));
        newpold = cell(1,length(ar.model(m).data));
        constVars = cell(1, length(ar.model(m).data));

        if(usePool)
            parfor d=1:length(ar.model(m).data)
                warning('off', 'symbolic:mupadmex:MuPADTextWarning');
                warning('off', 'symbolic:generate:FunctionNotVerifiedToBeValid')
                warning('off', 'symbolic:sym:sym:DeprecateExpressions')

                c = data(d).cLink;
                data_sym = arCalcData(config, model, data(d), m, c, d, doskip(d), matlab_version);
                newp{d} = data_sym.p;
                newfp{d} = data_sym.fp;
                newpold{d} = data_sym.pold;
                constVars{d} = data_sym.constVars;
                
                if(~doskip(d))
                    % header
                    fid_obsH = fopen([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.h'], 'W');
                    arWriteHFilesData(fid_obsH, data_sym);
                    fclose(fid_obsH);
                    movefile([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.h'], ...
                        [source_dir '/Compiled/' c_version_code '/' data(d).fkt '.h'],'f');
                    % body
                    fid_obs = fopen([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.c'], 'W');
                    arWriteCFilesData(fid_obs, matlab_version, config, m, c, d, data_sym);
                    fclose(fid_obs);
                    movefile([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.c'], ...
                        [source_dir '/Compiled/' c_version_code '/' data(d).fkt '.c'],'f');
                end
            end
        else
            for d=1:length(ar.model(m).data)
                c = data(d).cLink;
                try
                    data_sym = arCalcData(config, model, data(d), m, c, d, doskip(d), matlab_version);
                catch
                    data(d)
                end
                newp{d} = data_sym.p;
                newfp{d} = data_sym.fp;
                newpold{d} = data_sym.pold;
                constVars{d} = data_sym.constVars;
                
                if(~doskip(d))
                    % header
                    fid_obsH = fopen([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.h'], 'W');
                    arWriteHFilesData(fid_obsH, data_sym);
                    fclose(fid_obsH);
                    movefile([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.h'], ...
                        [source_dir '/Compiled/' c_version_code '/' data(d).fkt '.h'],'f');
                    % body
                    fid_obs = fopen([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.c'], 'W');
                    arWriteCFilesData(fid_obs, matlab_version, config, m, c, d, data_sym);
                    fclose(fid_obs);
                    movefile([source_dir '/Compiled/' c_version_code '/' data(d).fkt '_tmp.c'], ...
                        [source_dir '/Compiled/' c_version_code '/' data(d).fkt '.c'],'f');
                end
            end
        end

        % assign data
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).p = newp{d};
            ar.model(m).data(d).fp = newfp{d};
            ar.model(m).data(d).pold = newpold{d};
            ar.model(m).data(d).constVars = constVars{d};
        end
    else
        qdynparas = ismember(ar.model(m).p, ar.model(m).px) | ... %R2013a compatible
            ismember(ar.model(m).p, ar.model(m).pu); %R2013a compatible
        
        % conditions checksum
        checksum_cond = addToCheckSum(ar.model(m).fu);
        checksum_cond = addToCheckSum(ar.model(m).p(qdynparas), checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fv, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).N, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).cLink, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).z, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fz, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fp, checksum_cond);
        
        if isfield( ar.model(m), 'reducedForm' )
            % Store whether we express the model in totals or not
            for jt = 1 : numel( ar.model(m).removedStates )
                checksum_cond = addToCheckSum(ar.model(m).removedStates(jt).expressInTotal, checksum_cond);
            end
        end
        
        % append condition
        cindex = 1;
        
        ar.model(m).condition(cindex).status = 0;
        ar.model(m).condition(cindex).fu = ar.model(m).fu;
        ar.model(m).condition(cindex).fp = ar.model(m).fp(qdynparas);
        ar.model(m).condition(cindex).p = ar.model(m).p(qdynparas);
        ar.model(m).condition(cindex).checkstr = getCheckStr(checksum_cond);
        ar.model(m).condition(cindex).fkt = [ar.model(m).name '_' ar.model(m).condition(cindex).checkstr];
        ar.model(m).condition(cindex).dLink = [];
        
        % global checksum
        if(isempty(checksum_global))
            checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt);
        else
            checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt, checksum_global);
        end
        
        % skip calc conditions
        doskip = nan(1,length(ar.model(m).condition));
        for c=1:length(ar.model(m).condition)
            doskip(c) = ~forcedCompile && exist([source_dir '/Compiled/' ar.info.c_version_code '/' ar.model(m).condition(c).fkt '.c'],'file');
        end

        % calc conditions
        config = ar.config;
        model.name = ar.model(m).name;
        model.fv = ar.model(m).fv;
        model.fz = ar.model(m).fz;
        model.px0 = ar.model(m).px0;
        model.sym = ar.model(m).sym;
        model.t = ar.model(m).t;
        model.x = ar.model(m).x;
        model.u = ar.model(m).u;
        model.z = ar.model(m).z;
        model.us = ar.model(m).us;
        model.xs = ar.model(m).xs;
        model.vs = ar.model(m).vs;
        model.zs = ar.model(m).zs;
        model.N = ar.model(m).N;
        model.dvdx = ar.model(m).sym.dvdx;
        model.dvdu = ar.model(m).sym.dvdu;
        if isfield( ar.model(m), 'removedStates' )
            arFprintf( 2, 'Dealing with reduced model ...\n' );
            model.removedStates = ar.model(m).removedStates;
        end
        condition = ar.model(m).condition;
        newp = cell(1,length(ar.model(m).condition));
        newfp = cell(1,length(ar.model(m).condition));
        newpold = cell(1,length(ar.model(m).condition));
        newpx0 = cell(1,length(ar.model(m).condition));
        
        if(usePool)
            parfor c=1:length(ar.model(m).condition)         
                warning('off','symbolic:mupadmex:MuPADTextWarning');
                warning('off', 'symbolic:generate:FunctionNotVerifiedToBeValid')
                warning('off', 'symbolic:sym:sym:DeprecateExpressions')

                condition_sym = arCalcCondition(config, model, condition(c), m, c, doskip(c), matlab_version, simplifyEquations);
                newp{c} = condition_sym.p;
                newfp{c} = condition_sym.fp;
                newpold{c} = condition_sym.pold;
                newpx0{c} = condition_sym.px0;
                if(~doskip(c))
                    % header
                    fid_odeH = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], 'W'); % create header file
                    arWriteHFilesCondition(fid_odeH, config, condition_sym);
                    fclose(fid_odeH);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.h'],'f');
                    % body
                    fid_ode = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], 'W');
                    arWriteCFilesCondition(fid_ode, matlab_version, config, model, condition_sym, m, c, timedebug);
                    fclose(fid_ode);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.c'],'f');
                end
            end
        else        
            for c=1:length(ar.model(m).condition)
                condition_sym = arCalcCondition(config, model, condition(c), m, c, doskip(c), matlab_version, simplifyEquations);
                newp{c} = condition_sym.p;
                newfp{c} = condition_sym.fp;
                newpold{c} = condition_sym.pold;
                newpx0{c} = condition_sym.px0;
                if(~doskip(c))
                    % header
                    fid_odeH = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], 'W'); % create header file
                    arWriteHFilesCondition(fid_odeH, config, condition_sym);
                    fclose(fid_odeH);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.h'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.h'],'f');
                    % body
                    fid_ode = fopen([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], 'W');
                    arWriteCFilesCondition(fid_ode, matlab_version, config, model, condition_sym, m, c, timedebug);
                    fclose(fid_ode);
                    movefile([source_dir '/Compiled/' c_version_code '/' condition(c).fkt '_tmp.c'], ...
                        [source_dir '/Compiled/' c_version_code '/' condition(c).fkt '.c'],'f');
                end
            end
        end

        % assign conditions
        for c=1:length(ar.model(m).condition)
            ar.model(m).condition(c).p = newp{c};
            ar.model(m).condition(c).fp = newfp{c};
            ar.model(m).condition(c).pold = newpold{c};
            ar.model(m).condition(c).px0 = newpx0{c};
        end
        
        % plot setup
        if(~isfield(ar.model(m), 'plot'))
            ar.model(m).plot(1).name = ar.model(m).name;
        else
            ar.model(m).plot(end+1).name = ar.model(m).name;
        end
        ar.model(m).plot(end).doseresponse = false;
        ar.model(m).plot(end).dLink = 0;
        ar.model(m).plot(end).ny = 0;
        ar.model(m).plot(end).condition = {};
    end
end

ar.checkstr = getCheckStr(checksum_global);
ar.fkt = ['arSimuCalcFun_' ar.checkstr];

% write arSimuCalcFunctions
writeSimuCalcFunctions(debug_mode);

% compile
if ~forcedCompile && exist([pwd,filesep,ar.fkt,'.',mexext],'file')
    % do nothing
    fprintf('%s is already available: compiling skipped.\n',[ar.fkt,'.',mexext]);
elseif ( forcedCompile == 2 )
    arCompile(2, true, false, source_dir);
else
    arCompile(forcedCompile, false, false, source_dir);
end

% copy model and data backup files:
if isfield(ar.config,'backup_modelAndData') && ar.config.backup_modelAndData
    if ~isdir('./Models/Backup')
        mkdir('./Models/Backup');
    end
    if ~isdir(['./Models/Backup/',ar.checkstr])
        mkdir('./Models/Backup/',ar.checkstr);
    end
    if ~isdir('./Data/Backup')
        mkdir('./Data/Backup');
    end
    if ~isdir(['./Data/Backup/',ar.checkstr])
        mkdir('./Data/Backup/',ar.checkstr);
    end
    
    ar.setup.backup_data_folder = cell(size(ar.setup.datafiles));
    for i=1:length(ar.setup.datafiles)
        for j=1:length(ar.setup.datafiles{i})
            ar.setup.backup_data_folder{i}{j} = fullfile(pwd,['/Data/Backup/',ar.checkstr,'/']);% fullfile to prevent mixing of \ and /
            ar.setup.backup_data_folder_local{i}{j} = ['./Data/Backup/',ar.checkstr,'/'];
            [~,file,ext] = fileparts(ar.setup.datafiles{i}{j});
            source = ar.setup.datafiles{i}{j};
            target = [ar.setup.backup_data_folder{i}{j},file,ext];
            if ~isempty(source) && ~isempty(target) && strcmp(fullfile(strrep(source,pwd,'.')),fullfile(strrep(target,pwd,'.')))~=1
                copyfile(source,target);
            end
        end
    end
    ar.setup.backup_model_folder = cell(size(ar.setup.modelfiles));
    for i=1:length(ar.setup.modelfiles)
        if ~isempty(ar.setup.modelfiles{i})
            ar.setup.backup_model_folder{i} = fullfile(pwd,['/Models/Backup/',ar.checkstr,'/']);  % fullfile to prevent mixing of \ and /
            ar.setup.backup_model_folder_local{i} = ['./Models/Backup/',ar.checkstr,'/'];
            [~,file,ext] = fileparts(ar.setup.modelfiles{i});
            source = ar.setup.modelfiles{i};
            target = [ar.setup.backup_model_folder{i},file,ext];
            if ~isempty(source) && ~isempty(target) && strcmp(fullfile(strrep(source,pwd,'.')),fullfile(strrep(target,pwd,'.')))~=1
                copyfile(source,target);
            end
        end
    end
end


% link
arLink;

% refresh file cache
rehash

warning(warnreset);


% Calc Model
function arCalcModel(m, matlab_version, checkModel)
global ar

arFprintf(2, 'calculating model m%i, %s...\n', m, ar.model(m).name);

% make short strings
ar.model(m).xs = {};
ar.model(m).zs = {};
ar.model(m).us = {};
ar.model(m).vs = {};

for j=1:length(ar.model(m).x)
    ar.model(m).xs{j} = sprintf('x[%i]',j);
end
for j=1:length(ar.model(m).z)
    ar.model(m).zs{j} = sprintf('z[%i]',j);
end
for j=1:length(ar.model(m).u)
    ar.model(m).us{j} = sprintf('u[%i]',j);
end
for j=1:length(ar.model(m).fv)
    ar.model(m).vs{j} = sprintf('v[%i]',j);
end

% Obtain special functions (such as step functions)
specialFunc = ar.config.specialFunc;

% make syms
ar.model(m).sym.x = mySym(ar.model(m).x, specialFunc);
ar.model(m).sym.xs = mySym(ar.model(m).xs, specialFunc);
ar.model(m).sym.z = mySym(ar.model(m).z, specialFunc);
ar.model(m).sym.zs = mySym(ar.model(m).zs, specialFunc);
ar.model(m).sym.px0 = arMyStr2Sym(ar.model(m).px0);
ar.model(m).sym.u = mySym(ar.model(m).u, specialFunc);
ar.model(m).sym.us = mySym(ar.model(m).us, specialFunc);
ar.model(m).sym.v = mySym(strrep(strrep(ar.model(m).v, ' ', '_'),'.','_'), specialFunc);
ar.model(m).sym.vs = mySym(ar.model(m).vs, specialFunc);
ar.model(m).sym.fv = mySym(ar.model(m).fv, specialFunc);
[ ar.model(m).sym.dvdx, ar.model(m).sym.dvdu ] = cacheShortTerms( ar.model(m) );

% compartment volumes
if(~isempty(ar.model(m).pc)) 
    % make syms
    ar.model(m).sym.pc = arMyStr2Sym(ar.model(m).pc);
    ar.model(m).sym.C = arMyStr2Sym(ones(size(ar.model(m).N)));
    
    if(~isfield(ar.model(m),'isAmountBased') || ~ar.model(m).isAmountBased)
        for j=1:size(ar.model(m).N,1) % for every species j
            qinfluxwitheducts = ar.model(m).N(j,:) > 0 & sum(ar.model(m).N < 0,1) > 0;
            eductcompartment = zeros(size(qinfluxwitheducts));
            for jj=find(qinfluxwitheducts)
				eductcompartment(jj) = unique(ar.model(m).cLink(ar.model(m).N(:,jj)<0)); %R2013a compatible
            end
            
            cfaktor = arMyStr2Sym(ones(size(qinfluxwitheducts)));
            for jj=find(qinfluxwitheducts & eductcompartment~=ar.model(m).cLink(j))
                cfaktor(jj) = ar.model(m).sym.pc(eductcompartment(jj)) / ...
                    ar.model(m).sym.pc(ar.model(m).cLink(j));
            end
            ar.model(m).sym.C(j,:) = transpose(cfaktor);
        end
    else
        for j=1:size(ar.model(m).N,1) % for every species j
            ar.model(m).sym.C(j,:) = ar.model(m).sym.C(j,:) / ar.model(m).sym.pc(ar.model(m).cLink(j));
        end
    end
else
    ar.model(m).sym.C = arMyStr2Sym(ones(size(ar.model(m).N)));
end

% derivatives
if(~isempty(ar.model(m).sym.fv))
    ar.model(m).sym.dfvdx = myJacobian(ar.model(m).sym.fv, ar.model(m).sym.x);
    if(~isempty(ar.model(m).sym.us))
        ar.model(m).sym.dfvdu = myJacobian(ar.model(m).sym.fv, ar.model(m).sym.u);
    else
        ar.model(m).sym.dfvdu = arMyStr2Sym(ones(length(ar.model(m).sym.fv), 0));
    end
else
    ar.model(m).sym.dfvdx = arMyStr2Sym(ones(0, length(ar.model(m).sym.x)));
    ar.model(m).sym.dfvdu = arMyStr2Sym(ones(0, length(ar.model(m).sym.u)));
end

if length(ar.model(m).sym.x)<100 && length(ar.model(m).p)<500
    ar.model(m).qdvdx_nonzero = logical(ar.model(m).sym.dfvdx~=0);
    ar.model(m).qdvdu_nonzero = logical(ar.model(m).sym.dfvdu~=0);
    
    if ( checkModel )
        tmpsym = ar.model(m).sym.dfvdx;
        tmpsym = arSubs(tmpsym, ar.model(m).sym.x,  rand(size(ar.model(m).sym.x)), matlab_version);
        tmpsym = arSubs(tmpsym, ar.model(m).sym.u,  rand(size(ar.model(m).sym.u)), matlab_version);
        tmpsym = arSubs(tmpsym, arMyStr2Sym(ar.model(m).p), rand(size(ar.model(m).p)), matlab_version);

        try
            ar.model(m).q_negative = double(tmpsym) < 0;
        catch ERR
            for i=1:length(tmpsym(:))
                try
                    double(tmpsym(i));
                catch
                    warning('Failed to check dfvdx due to use of symbolic function:')
                    tmpsym(i)
                end
            end
        end
    end
    
    tmpsym = ar.model(m).sym.dfvdu;
    tmpsym = arSubs(tmpsym, ar.model(m).sym.x, rand(size(ar.model(m).sym.x)), matlab_version);
    tmpsym = arSubs(tmpsym, ar.model(m).sym.u, rand(size(ar.model(m).sym.u)), matlab_version);
    tmpsym = arSubs(tmpsym, arMyStr2Sym(ar.model(m).p), rand(size(ar.model(m).p)), matlab_version);
    
    ar.model(m).qdvdu_negative = double(tmpsym) < 0;
    
    
    tmpzeros = (ar.model(m).N .* ar.model(m).sym.C) * ar.model(m).sym.dfvdx;
    ar.model(m).nnz = arnansum(arnansum(logical(tmpzeros~=0))) + arnansum(arnansum(logical(tmpzeros~=0))==0);
    if(length(ar.model(m).x) * log(length(ar.model(m).x)) > ar.model(m).nnz)
        ar.config.useSparseJac = 1;
    end
else
    arFprintf(2, 'Model m%i too large (%i variables, %i parameters)... skipping check of model structure\n', m, length(ar.model(m).sym.x), length(ar.model(m).p));
    ar.config.useSparseJac = 1;
    tmpzeros = (ar.model(m).N .* ar.model(m).sym.C) * ar.model(m).sym.dfvdx;
    ar.model(m).nnz = arnansum(arnansum(logical(tmpzeros~=0))) + arnansum(arnansum(logical(tmpzeros~=0))==0);
end



% Calc Condition
function condition = arCalcCondition(config, model, condition, m, c, doskip, matlab_version, simplifyEquations)

if(doskip)
    arFprintf(2, 'calculating condition m%i c%i, %s...skipped\n', m, c, model.name);
else
    arFprintf(2, 'calculating condition m%i c%i, %s...\n', m, c, model.name);
end

% hard code conditions
specialFunc = config.specialFunc;
condition.sym.p = arMyStr2Sym(condition.p);
try
    condition.sym.fp = arMyStr2Sym(condition.fp);
catch
    error( invalidSym( 'Invalid expression in condition transformations', condition.fp, condition.p ) );
end

% Were model reductions applied? Make sure to transfer the initial conditions correctly
if isfield( model, 'removedStates' )
    condition = addPoolSubstitutions( condition, model.removedStates );
end

condition.sym.fpx0 = arMyStr2Sym(model.px0);
condition.sym.fpx0 = arSubs(condition.sym.fpx0, condition.sym.p, condition.sym.fp, matlab_version);
condition.sym.fv = mySym(model.fv, specialFunc);
condition.sym.fv = arSubs(condition.sym.fv, condition.sym.p, condition.sym.fp, matlab_version);
condition.sym.fu = mySym(condition.fu, specialFunc);
condition.sym.fu = arSubs(condition.sym.fu, condition.sym.p, condition.sym.fp, matlab_version);
condition.sym.fz = mySym(model.fz, specialFunc);
condition.sym.fz = mysubsrepeated(condition.sym.fz, model.sym.z, condition.sym.fz, matlab_version); % Substitute references to derived variables
condition.sym.fz = arSubs(condition.sym.fz, condition.sym.p, condition.sym.fp, matlab_version);
condition.sym.C = arSubs(model.sym.C, condition.sym.p, condition.sym.fp, matlab_version);

% Replace inline arrays (e.g. [3,4,2,5,2] with variable names, declaring
% the array elsewhere as a static const (used for the fixed input spline)
constVars = {};
cVars = {};
for i=1:length(condition.sym.fu)
    [ tmp, constVars, cVars ] = replaceFixedArrays( char(condition.sym.fu(i)), constVars, cVars, c, 'condi' );
    condition.sym.fu(i) = arMyStr2Sym( tmp );
end


for i=1:length(condition.sym.fv)
    [ tmp, constVars, cVars ] = replaceFixedArrays( char(condition.sym.fv(i)), constVars, cVars, c, 'condi' );
    condition.sym.fv(i) = arMyStr2Sym( tmp );
end

condition.constVars = sprintf( '%s;\n', constVars{:} );

% predictor
condition.sym.fv = arSubs(condition.sym.fv, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);
condition.sym.fu = arSubs(condition.sym.fu, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);
condition.sym.fz = arSubs(condition.sym.fz, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);

% remaining initial conditions
varlist = symvar(condition.sym.fpx0);
condition.px0 = sym2str(varlist);

% remaining parameters
varlist = union( symvar([condition.sym.fv(:); condition.sym.fu(:); condition.sym.fz(:); condition.sym.fpx0(:)]), symvar( condition.sym.C ) );
condition.pold = condition.p;
condition.p = setdiff(setdiff(setdiff(setdiff(setdiff(sym2str(varlist), model.x), model.u), model.z), 't'), cVars);
condition.dfxdx_rowVals = [];
condition.dfxdx_colptrs = [];  

% We cannot skip compilation when we're using the faster splines. The
% reason is that it is impossible to determine how many splines end up in
% the derivatives without determining the derivatives.
if ( config.turboSplines == 1 )
    doskip = false;
end

if(doskip)
    condition.ps = {};
    condition.qfu_nonzero = [];
    condition.qdvdx_nonzero = [];
    condition.qdvdu_nonzero = [];
    condition.qdvdp_nonzero = [];
    condition.dvdx = {};
    condition.dvdu = {};
    condition.dvdp = {};
    condition.qdfxdx_nonzero = [];     
    condition.dfxdx = {};
    condition.su = {};
    condition.sx = {};
    condition.qfsv_nonzero = [];
    condition.sv = {};
    condition.sz = {};    
    condition.splines = [];
    condition.constVars = '';
    
    return;
end

if (simplifyEquations > 0)
    condition.sym.fz = simplify( condition.sym.fz );
end
if (simplifyEquations > 1)
    condition.sym.fv = simplify( condition.sym.fv );
end


% make short strings
condition.ps = {};
for j=1:length(condition.p)
    condition.ps{j} = sprintf('p[%i]',j);
end

% make syms
condition.sym.p = arMyStr2Sym(condition.p);
condition.sym.ps = arMyStr2Sym(condition.ps);
condition.sym.px0s = arSubs(arMyStr2Sym(condition.px0), ...
    condition.sym.p, condition.sym.ps, matlab_version);

% make syms
condition.sym.fv = arSubs(condition.sym.fv, model.sym.x, model.sym.xs, matlab_version);
condition.sym.fv = arSubs(condition.sym.fv, model.sym.u, model.sym.us, matlab_version);
condition.sym.fv = arSubs(condition.sym.fv, condition.sym.p, condition.sym.ps, matlab_version);

condition.sym.fu = arSubs(condition.sym.fu, condition.sym.p, condition.sym.ps, matlab_version);

condition.sym.fz = arSubs(condition.sym.fz, model.sym.x, model.sym.xs, matlab_version);
condition.sym.fz = arSubs(condition.sym.fz, model.sym.u, model.sym.us, matlab_version);
condition.sym.fz = arSubs(condition.sym.fz, condition.sym.p, condition.sym.ps, matlab_version);

condition.sym.fpx0 = arSubs(condition.sym.fpx0, condition.sym.p, condition.sym.ps, matlab_version);

% remove zero inputs
condition.qfu_nonzero = logical(condition.sym.fu ~= 0);
if(~isempty(model.sym.us))
    condition.sym.fv = arSubs(condition.sym.fv, model.sym.us(~condition.qfu_nonzero), ...
        arMyStr2Sym(zeros(1,sum(~condition.qfu_nonzero))), matlab_version);
end

% derivatives
if(~isempty(condition.sym.fv))
    condition.sym.dfvdx = myJacobian(condition.sym.fv, model.sym.xs);
    if(~isempty(model.sym.us))
        condition.sym.dfvdu = myJacobian(condition.sym.fv, model.sym.us);
    else
        condition.sym.dfvdu = arMyStr2Sym(ones(length(condition.sym.fv), 0));
    end
    condition.sym.dfvdp = myJacobian(condition.sym.fv, condition.sym.ps);
else
    condition.sym.dfvdx = arMyStr2Sym(ones(0, length(model.sym.xs)));
    condition.sym.dfvdu = arMyStr2Sym(ones(0, length(model.sym.us)));
    condition.sym.dfvdp = arMyStr2Sym(ones(0, length(condition.sym.ps)));
end

% flux signs
condition.qdvdx_nonzero = logical(condition.sym.dfvdx~=0);
condition.qdvdu_nonzero = logical(condition.sym.dfvdu~=0);
condition.qdvdp_nonzero = logical(condition.sym.dfvdp~=0);

% short terms
condition.sym.dvdx = model.dvdx;
if ~isempty(condition.qdvdx_nonzero)
    condition.sym.dvdx(~condition.qdvdx_nonzero) = 0;
end
condition.sym.dvdu = model.dvdu;
if ~isempty(condition.qdvdu_nonzero)
    condition.sym.dvdu(~condition.qdvdu_nonzero) = 0;
end
    
condition.sym.dvdp = arMyStr2Sym(zeros(length(model.vs), length(condition.ps)));
for j=1:length(model.vs)
    for i=1:length(condition.ps)
        if(condition.qdvdp_nonzero(j,i))
            condition.sym.dvdp(j,i) = arMyStr2Sym(sprintf('dvdp[%i]', j + (i-1)*length(model.vs)));
        end
    end
end

% do we have variable volumes?
if ( ~isempty( symvar( condition.sym.C ) ) )
    warning( 'Variable volume detected. Sensitivities for variable volumes are still in an experimental stage' );
    
    for a = 1 : length( condition.sym.p )
        condition.sym.dfcdp(:,a) = (diff(model.N.*condition.sym.C, condition.sym.p(a)))*condition.sym.fv;
    end
    condition.sym.dfcdp = arSubs(condition.sym.dfcdp, condition.sym.p, condition.sym.ps, matlab_version);
end

% make equations
condition.sym.C = arSubs(condition.sym.C, condition.sym.p, condition.sym.ps, matlab_version);
condition.sym.fx = (model.N .* condition.sym.C) * transpose(model.sym.vs);
firstcol = true;
% Jacobian dfxdx
if(config.useSensis || config.useJacobian)
    condition.sym.dfxdx = (model.N .* condition.sym.C) * condition.sym.dvdx;
    condition.qdfxdx_nonzero = logical(condition.sym.dfxdx~=0);
    condition.sym.dfxdx_nonzero = arMyStr2Sym(zeros(1, arnansum(arnansum(condition.qdfxdx_nonzero))));
    for j=1:length(model.xs)
        for i=1:length(model.xs)
            if(i==1)
               firstcol = true;
            end
            if(~isempty(condition.qdfxdx_nonzero) && condition.qdfxdx_nonzero(i,j))
                condition.dfxdx{i,j} = sprintf('dfxdx[%i]', i + (j-1)*length(model.xs));                
                condition.dfxdx_rowVals = [condition.dfxdx_rowVals i-1];
                condition.sym.dfxdx_nonzero(length(condition.dfxdx_rowVals)) = condition.sym.dfxdx(i,j);
                if(firstcol)
                    condition.dfxdx_colptrs = [condition.dfxdx_colptrs length(condition.dfxdx_rowVals)-1];
                    firstcol = false;
                end                
            else
                condition.dfxdx{i,j} = '0';
            end
            if(firstcol && i==length(model.xs) && ( size(condition.sym.dfxdx_nonzero, 2) > 0) )
                condition.dfxdx_rowVals = [condition.dfxdx_rowVals i-1];
                condition.sym.dfxdx_nonzero(length(condition.dfxdx_rowVals)) = arMyStr2Sym('RCONST(0.0)');
                condition.dfxdx_colptrs = [condition.dfxdx_colptrs length(condition.dfxdx_rowVals)-1];
                firstcol = false;
            end
        end
    end
end
condition.dfxdx_colptrs = [condition.dfxdx_colptrs length(condition.dfxdx_rowVals)];
condition.sym.dfzdu = myJacobian(condition.sym.fz, model.sym.us);
condition.sym.dfzdx = myJacobian(condition.sym.fz, model.sym.xs);

% sensitivities
if(config.useSensis)
	% su
    condition.su = cell(length(model.us), 1);
    for j=1:length(model.us)
        if(condition.qfu_nonzero(j))
            condition.su{j} = sprintf('su[%i]', j);
        else
            condition.su{j} = '0';
        end
    end
    condition.sym.su = mySym(condition.su, specialFunc);
    
    % input derivatives 
    if(~isempty(condition.sym.ps))
        if(~isempty(condition.sym.fu))
            condition.sym.dfudp = ...
                myJacobian(condition.sym.fu, condition.sym.ps);
        else
            condition.sym.dfudp = arMyStr2Sym(ones(0,length(condition.sym.ps)));
        end
        
        % Replace spline derivatives
        for j = 1 : length( model.u )
            if ( strfind( condition.fu{j}, 'spline' ) )
                for j2 = 1 : length( condition.sym.dfudp(j,:) )
                    ustr = char(condition.sym.dfudp(j, j2));
                    ustr = repSplineDer( ustr );
                    condition.sym.dfudp(j,j2) = arMyStr2Sym(ustr);
                end
            end
        end
                
        % This function checks whether the inputs were sensible and
        % gives a warning for problematic discontinuities in the sensitivities.
        % for j = 1 : length( model.u )
        %     condition.sym.dfudp(j,:) = verifyRow( condition.sym.dfudp(j,:), condition.sym.fu(j), 'input' );
        % end
    end
    
	% sx
    condition.sx = cell(length(model.xs), 1);
    for j=1:length(model.xs)
        condition.sx{j} = sprintf('sx[%i]', j);
    end
	condition.sym.sx = arMyStr2Sym(condition.sx);
    
    condition.sym.fsv1 = condition.sym.dvdx * condition.sym.sx + ...
        condition.sym.dvdu * condition.sym.su;
    % fsv2 = condition.sym.dvdp;
    
	% sv
    condition.sv = cell(length(model.vs), 1);
    for j=1:length(model.vs)
        condition.sv{j} = sprintf('sv[%i]', j);
    end
    condition.sym.sv = arMyStr2Sym(condition.sv);
    
    if(config.useSensiRHS)
        condition.sym.fsx = (model.N .* condition.sym.C) * condition.sym.sv;
    end
    
    % sx initials
    if(~isempty(condition.sym.fpx0))
        condition.sym.fsx0 = myJacobian(condition.sym.fpx0, condition.sym.ps);
    else
        condition.sym.fsx0 = arMyStr2Sym(ones(0, length(condition.sym.ps)));
    end
    
    % partial derivative right hand side w.r.t. parameters (keeping x fixed)
    if(isfield(condition.sym, 'dfudp'))
        condition.sym.dfxdp = (model.N .* condition.sym.C) * (condition.sym.dvdp + condition.sym.dvdu * condition.sym.dfudp);
    else
        condition.sym.dfxdp = (model.N .* condition.sym.C) * (condition.sym.dvdp);
    end
    
    % initial sensitivities used for steady state constraints (total derivative at initial condition)
    if(isfield(condition.sym, 'dfudp'))
        condition.sym.dfxdp0 = (model.N .* condition.sym.C) * (condition.sym.dvdp + ...
            condition.sym.dvdx*condition.sym.fsx0 + ...
            condition.sym.dvdu * condition.sym.dfudp);
    else
        condition.sym.dfxdp0 = (model.N .* condition.sym.C) * (condition.sym.dvdp + ...
            condition.sym.dvdx*condition.sym.fsx0);
    end

    % Add variable volume terms here
    if (isfield(condition.sym, 'dfcdp'))
        condition.sym.dfxdp0 = condition.sym.dfxdp0 + condition.sym.dfcdp;
        condition.sym.dfxdp  = condition.sym.dfxdp  + condition.sym.dfcdp;
    end
    
    % derivatives fz
    condition.sym.dfzdp = myJacobian(condition.sym.fz, condition.sym.ps);
    
    % sz
    condition.sz = cell(length(model.zs), 1);
    for j=1:length(model.zs)
        condition.sz{j,1} = sprintf('sz[%i]', j);
    end
    condition.sym.sz = arMyStr2Sym(condition.sz);
    
    condition.sym.fsz1 = condition.sym.dfzdx * condition.sym.sx;
    if(isfield(condition.sym, 'dfudp'))
        condition.sym.fsz2 = condition.sym.dfzdu * condition.sym.dfudp + ...
            condition.sym.dfzdp;
    else
        condition.sym.fsz2 = condition.sym.dfzdp;
    end
end

% Also generate second order equations (for second order steady state
% sensitivities)
%condition.sym.d2fxdxdp = jacobian(condition.sym.dfxdx, condition.p);
%condition.sym.d2fxdpdq = jacobian(condition.sym.dfxdp, condition.p);

% We need room for the spline coefficients if we use turbosplines
if ( config.turboSplines == 1 )
    condition = uniqueSplines( condition );
else
    condition.splines = [];
end

% Replaces arrays with static variables
% And outputs a list of the new static variables created.
function [ finalStr, static_variables, cvars ] = replaceFixedArrays( str, static_variables, cvars, condiID, name )

    idx = numel(cvars);
    finalStr = str;
    [strOut, values] = findArrays( str );
        
    % Generate static variables for all the const arrays
    for a = 1 : numel( strOut )
        idx = idx + 1;
        vString = sprintf('%.25e, ', values{a});
        vString = vString(1:end-2);
        cvars{end+1} = sprintf( 'static_chunk_%s_%d_%d', name, condiID, idx ); %#ok
        static_variables{end+1} = sprintf( 'static const double static_chunk_%s_%d_%d[%d] = { %s }', name, condiID, idx, numel( values{a} ), vString ); %#ok
        finalStr = strrep( finalStr, strOut{a}, sprintf( 'static_chunk_%s_%d_%d', name, condiID, idx ) );
    end

% Grab arrays from string ([1,2,3,...5,6,7])
function [strOut, values] = findArrays( str )
    %loc1 = regexp( str, '\[([ ()0123456789\.\,\^\*\-\+\/]+)[,]([ ()0123456789\.\,\^\*\-\+\/]+)\]' );
    %loc2 = regexp( str, '([ ()0123456789\.\,\^\*\-\+\/]+)[,]([ ()0123456789\.\,\^\*\-\+\/]+)\]', 'end' );
    
    % This is less robust than the regular expression based system, but
    % orders of magnitude faster.
    loc1 = str=='[';
    loc2 = str==']';
       
    % Encode brackets so we can find closed pairs
    empty = zeros(numel(str), 1);
    empty(loc1) =  1;
    empty(loc2) = -1;
    nz = find(abs(empty)>0);
    brackets = empty(nz);
    dbrackets = diff(brackets);
    % Brackets with a derivative of -2 are closed lowest level brackets
    L = find( dbrackets == -2 );
    loc1 = nz(L);
    loc2 = nz(L+1);
    
    a = 1;
    strOut = cell( 1, numel(loc1) );
    values = cell( 1, numel(loc1) );
    for i = 1 : numel( loc1 )
        curValues = arrayVal( str( loc1(i) + 1 : loc2(i) - 1 ) );
        if ~isempty( curValues )
            strOut{a} = str( loc1(i) : loc2(i) );
            values{a} = curValues;
            a = a + 1;
        end
    end
    strOut = strOut(1:a-1);
    values = values(1:a-1);
    

% Grab values from flat array and force them to doubles
function vals = arrayVal( str )
    spl = strsplit( str, ',' );
    try
        vals = cellfun( @str2num, spl );
        
        % Not actually an array! Do not replace!
        if ( numel( vals ) == 1 )
            vals = [];
        end
    catch
        vals = [];
        %error( 'Found a malformed array: [%s]. Only numeric values are allowed!', str )
    end

% This function checks whether any deltas appear in the sensitivity
% equations. They are incompatible with continuous optimizers
function sensBlock = verifyRow( sensBlock, func, location )
    for k = 1 : size( sensBlock, 2 )
        jacElemStr = char(sensBlock(1,k));
        if (numel(strfind(jacElemStr, 'dirac(')>0))
            % Attempt to simplify (sometimes this gets rid of dirac
            % functions
            fprintf( 'Dirac detected in derivative => attempting to simplify ...' );
            try
                jacElem = simplify( sensBlock(1,k) );
                jacElemStr = char(jacElem);
            catch
                pretty( sensBlock(1,k) );
                error( 'Sensitivity expression for %s equation contains a non-compilable error (e.g. division by zero)', char(location) );
            end

            if (numel(strfind(jacElemStr, 'dirac(')>0))
                fprintf( ' [FAILED]\n' );
                % We failed to resolve the derivatives. Give up, but
                % let the user know the offending lines.
                message = {     'UNRESOLVABLE DERIVATIVE FOUND IN SENSITIVITY JACOBIAN\n\n'         , ...
                                'Equation (in ', char(location), ' section):\n\n', char( func )     , ...
                                '\n\nSensitivity equation:\n\n', char( sensBlock(1,k) )             , ...
                                '\n\nThis is likely due to a step function with variable location'  , ...
                                '\nparameter. Consider changing the equation to a continuous\n'     , ...
                                'function (e.g. smoothstep1).\n\nSetting corresponding sensitivity ' , ...
                                'to zero.\n\nThis means the sensitivity solution is now incorrect. ', ...
                                'Which\nmeans the corresponding parameter cannot be optimized using', ...
                                ' a\nderivative based optimization algorithm.\n\n'                   , ...
                                'Hit any key to proceed compiling (at your own risk).'};
                warning( sprintf( sprintf( '%s', message{:} ) ) ); %#ok
            else
                % This fixed it. Substitute Jacobian entry.
                fprintf( ' [ OK ]\n' );
                sensBlock(1,k) = jacElem;
            end
        end
    end            
    
% Calc Data
function data = arCalcData(config, model, data, m, c, d, doskip, matlab_version)

if(doskip)
    arFprintf(2, 'calculating data m%i d%i -> c%i, %s...skipped\n', m, d, c, data.name);
else
    arFprintf(2, 'calculating data m%i d%i -> c%i, %s...\n', m, d, c, data.name);
end

% Grab special functions list
specialFunc = config.specialFunc;

% hard code conditions
data.sym.p  = arMyStr2Sym(data.p);
data.sym.fp = arMyStr2Sym(data.fp);

% Were model reductions applied? Make sure to transfer the initial conditions correctly
if isfield( model, 'removedStates' )
    data = addPoolSubstitutions( data, model.removedStates );
end

% Replace inline arrays (e.g. [3,4,2,5,2] with variable names, declaring
% the array elsewhere as a static const (used for the fixed input spline)
constVars = {}; cVars = [];
for jy = 1 : numel( data.fy )
    [ data.fy{jy}, constVars, cVars ] = replaceFixedArrays( data.fy{jy}, constVars, {}, c, 'data' );
end
data.constVars = sprintf( '%s;\n', constVars{:} );

data.sym.fy = mySym(data.fy, specialFunc);
data.sym.fy = arSubs(data.sym.fy, model.sym.v, model.sym.fv, matlab_version);
data.sym.fy = arSubs(data.sym.fy, data.sym.p, data.sym.fp, matlab_version);
try
    data.sym.fystd = arMyStr2Sym(data.fystd);
catch
    error( 'Invalid expression in error model\n%s', [sprintf('%s', data.fystd{:} )] ); %#ok
end
data.sym.fystd = arSubs(data.sym.fystd, model.sym.v, model.sym.fv, matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, data.sym.p, data.sym.fp, matlab_version);

data.sym.fu = mySym(data.fu, specialFunc);
data.sym.fu = arSubs(data.sym.fu, data.sym.p, data.sym.fp, matlab_version);
data.qfu_nonzero = logical(data.sym.fu ~= 0);

% predictor
data.sym.fu = arSubs(data.sym.fu, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);
data.sym.fy = arSubs(data.sym.fy, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);

data.sym.fy = arSubs(data.sym.fy, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, arMyStr2Sym(model.t), arMyStr2Sym('t'), matlab_version);

% remaining parameters
varlist = symvar([data.sym.fy(:); data.sym.fystd(:)]);
data.pold = data.p;
othervars = union(union(union(union(model.x, model.u), model.z), data.y), 't');
data.p = setdiff( setdiff(union(sym2str(varlist), data.p_condition), othervars), cVars); %R2013a compatible

% Union's behaviour is different when first arg is empty. In this case, a
% flip of the parameter vector is typically required.
if ( size( data.p, 1 ) ~= 1 )
    data.p = data.p.';
end

if(doskip)
    data.ps = {};
    data.ys = [];
    data.qu_measured = [];
    data.qx_measured = [];
    data.dfydxnon0 = {};
    data.sx = {};
    data.sz = {};
    data.sy = {};
    
    return;
end

% make short strings
data.ps = cell(size(data.p));
for j=1:length(data.p)
    data.ps{j} = sprintf('p[%i]',j);
end
data.ys = {};
for j=1:length(data.y)
    data.ys{j} = sprintf('y[%i]',j);
end

% make syms
data.sym.p = arMyStr2Sym(data.p);
data.sym.ps = arMyStr2Sym(data.ps);
data.sym.y = arMyStr2Sym(data.y);
data.sym.ys = arMyStr2Sym(data.ys);

% substitute
data.sym.fy = arSubs(data.sym.fy, ...
    model.sym.x, model.sym.xs, matlab_version);
data.sym.fy = arSubs(data.sym.fy, ...
    model.sym.u, model.sym.us, matlab_version);
data.sym.fy = arSubs(data.sym.fy, ...
    model.sym.z, model.sym.zs, matlab_version);
data.sym.fy = arSubs(data.sym.fy, ...
    data.sym.p, data.sym.ps, matlab_version);

data.sym.fystd = arSubs(data.sym.fystd, ...
    model.sym.x, model.sym.xs, matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, ...
    model.sym.u, model.sym.us, matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, ...
    model.sym.z, model.sym.zs, matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, ...
    data.sym.y, data.sym.ys, matlab_version);
data.sym.fystd = arSubs(data.sym.fystd, ...
    data.sym.p, data.sym.ps, matlab_version);

% derivatives fy
if(~isempty(data.sym.fy))
    if(~isempty(model.sym.us))
        data.sym.dfydu = myJacobian(data.sym.fy, model.sym.us);
    else
        data.sym.dfydu = arMyStr2Sym(ones(length(data.y), 0));
    end
    if(~isempty(model.x))
        data.sym.dfydx = myJacobian(data.sym.fy, model.sym.xs);
    else
        data.sym.dfydx = [];
    end
    if(~isempty(model.z))
        data.sym.dfydz = myJacobian(data.sym.fy, model.sym.zs);
    else
        data.sym.dfydz = [];
    end
	data.sym.dfydp = myJacobian(data.sym.fy, data.sym.ps);
else
	data.sym.dfydu = [];
	data.sym.dfydx = [];
    data.sym.dfydz = [];
	data.sym.dfydp = [];
end

% what is measured ?
data.qu_measured = sum(logical(data.sym.dfydu~=0),1)>0;
data.qx_measured = sum(logical(data.sym.dfydx~=0),1)>0;
data.qz_measured = sum(logical(data.sym.dfydz~=0),1)>0;

% derivatives fystd
if(~isempty(data.sym.fystd))
    if(~isempty(model.sym.us))
        data.sym.dfystddu = myJacobian(data.sym.fystd, model.sym.us);
    else
        data.sym.dfystddu = arMyStr2Sym(ones(length(data.y), 0));
    end
    if(~isempty(model.x))
        data.sym.dfystddx = myJacobian(data.sym.fystd, model.sym.xs);
    else
        data.sym.dfystddx = [];
    end
    if(~isempty(model.z))
        data.sym.dfystddz = myJacobian(data.sym.fystd, model.sym.zs);
    else
        data.sym.dfystddz = [];
    end
    data.sym.dfystddp = myJacobian(data.sym.fystd, data.sym.ps);
    data.sym.dfystddy = myJacobian(data.sym.fystd, data.sym.ys);
else
    data.sym.dfystddu = [];
    data.sym.dfystddp = [];
	data.sym.dfystddy = [];
    data.sym.dfystddx = [];
    data.sym.dfystddz = [];
end

% observed directly and exclusively
data.dfydunon0 = logical(data.sym.dfydu ~= 0);
data.dfydxnon0 = logical(data.sym.dfydx ~= 0);
data.dfydznon0 = logical(data.sym.dfydz ~= 0);

% dzdx sensitivities
data.dfzdx = cell(length(model.zs), length(model.xs));
for j=1:length(model.zs)
    for i=1:length(model.xs)
        data.dfzdx{j,i} = sprintf('dfzdx[%i]', j + (i-1)*length(model.zs));
    end
end
data.sym.dfzdx = arMyStr2Sym(data.dfzdx);

% calculate y_scale
if(isempty(data.sym.dfydx) && isempty(data.sym.dfydz))
    data.sym.y_scale = arMyStr2Sym(zeros(length(data.sym.fy)));
else
    if(~isempty(data.sym.dfydx))
        data.sym.y_scale = data.sym.dfydx;    
    end
    if(~isempty(data.sym.dfydz))
        tmpfsz = data.sym.dfydz * ...
            data.sym.dfzdx;           
        if(~isempty(data.sym.dfydx))
            data.sym.y_scale = data.sym.y_scale + tmpfsz;
        else   
            data.sym.y_scale = tmpfsz;
        end
    end       
end

if(config.useSensis)
    % sx sensitivities
    data.sx = {};
    for j=1:length(model.xs)
        for i=1:length(data.p_condition)
            data.sx{j,i} = sprintf('sx[%i]', j + (i-1)*length(model.xs));
        end
    end
    data.sym.sx = arMyStr2Sym(data.sx);
    
    % su
    data.su = cell(length(model.us), length(data.p_condition));
    for j=1:length(model.us)
        for i=1:length(data.p_condition)
            if(data.qfu_nonzero(j))
                data.su{j,i} = sprintf('su[%i]', j + (i-1)*length(model.us));
            else
                data.su{j,i} = '0';
            end
        end
    end
    data.sym.su = arMyStr2Sym(data.su);
    
    % sz
    data.sz = cell(length(model.zs), length(data.p_condition));
    for j=1:length(model.zs)
        for i=1:length(data.p_condition)
            data.sz{j,i} = sprintf('sz[%i]', j + (i-1)*length(model.zs));
        end
    end
    data.sym.sz = arMyStr2Sym(data.sz);
    
    % sy sensitivities
    data.sy = {};
    for j=1:length(data.sym.fy)
        for i=1:length(data.sym.ps)
            data.sy{j,i} = sprintf('sy[%i]', j + (i-1)*length(data.sym.fy));
        end
    end
	data.sym.sy = arMyStr2Sym(data.sy);
    
    % parameters that appear in conditions
    if(~isempty(data.p_condition))
        qdynpara = ismember(data.p, data.p_condition); %R2013a compatible
    else
        qdynpara = false(size(data.p));
    end
    
    % calculate sy
    if(~isempty(data.sym.sy))
        data.sym.fsy = data.sym.dfydp;
        
        if(~isempty(data.p_condition))
            tmpfsx = data.sym.dfydx * ...
                data.sym.sx;
            tmpfsu = data.sym.dfydu * ...
                data.sym.su;
            tmpfsz = data.sym.dfydz * ...
                data.sym.sz;
            if(~isempty(model.u))
                data.sym.fsy(:,qdynpara) = data.sym.fsy(:,qdynpara) + tmpfsu;
            end
            if(~isempty(model.x))
                data.sym.fsy(:,qdynpara) = data.sym.fsy(:,qdynpara) + tmpfsx;
            end
            if(~isempty(model.z))
                data.sym.fsy(:,qdynpara) = data.sym.fsy(:,qdynpara) + tmpfsz;
            end
        end
    else
        data.sym.fsy = [];
    end
    
    % calculate systd
    if(~isempty(data.sym.sy))
        data.sym.fsystd = data.sym.dfystddp + ...
            data.sym.dfystddy * data.sym.sy;
                
        if(~isempty(data.p_condition))
            tmpfsx = data.sym.dfystddx * ...
                data.sym.sx;
            tmpfsu = data.sym.dfystddu * ...
                data.sym.su;
            tmpfsz = data.sym.dfystddz * ...
                data.sym.sz;
            if(~isempty(model.u))
                data.sym.fsystd(:,qdynpara) = data.sym.fsystd(:,qdynpara) + tmpfsu;
            end
            if(~isempty(model.x))
                data.sym.fsystd(:,qdynpara) = data.sym.fsystd(:,qdynpara) + tmpfsx;
            end
            if(~isempty(model.z))
                data.sym.fsystd(:,qdynpara) = data.sym.fsystd(:,qdynpara) + tmpfsz;
            end
        end
    else
        data.sym.fsystd = [];
    end
end

function condition = addPoolSubstitutions( condition, removedStates )

    % Substitute current state initials into total pool
    for jp = 1 : numel( removedStates )
        totalPool = subs( removedStates(jp).totalPoolStates, condition.sym.p.', condition.sym.fp );
        totalPool_str = char( totalPool );

        % Remove the original initial condition parameters
        condition.fp( strcmp( condition.p, removedStates(jp).initial ) ) = [];
        condition.p(  strcmp( condition.p, removedStates(jp).initial ) ) = [];

        % Remove the original initial condition parameters
        condition.sym.fp( find( condition.sym.p == removedStates(jp).sym.initial ) ) = []; %#ok
        condition.sym.p(  find( condition.sym.p == removedStates(jp).sym.initial ) ) = []; %#ok

        if ( removedStates(jp).expressInTotal == 0 )
            condition.fp{ strcmp( condition.p, removedStates(jp).totalVariable ) }             = totalPool_str;
            condition.sym.fp( find( condition.sym.p == removedStates(jp).sym.totalVariable ) ) = totalPool; %#ok
        end
    end

% substitute until no more changes (for self-substitutions of derived
% variables)
function out = mysubsrepeated(in, old, new, matlab_version)
    done = false;
    
    k = 0; orig = in;
    while ( ~done )
        out = arSubs(in, old, new, matlab_version);
        
        if ( k > 5 )
            v = '';
            for c = 1 : length( orig )
                if ~isequal( in(c), out(c) )
                    v = sprintf( '%s\n%s = %s', v, char(old(c)), char(orig(c)) );
                end
            end
            s = sprintf( 'Substitution recursion limit (5) exceeded!\nSolutions that cannot be obtained by simple substitution are not supported.\nDo you have any cyclic substitutions?\n%s\n', v );
            
            error( s ); %#ok
        end        
        
        % No more changes?
        if ( isequal( out, in ) == 1 )
            done = true;
        else
            in = out;
        end
        k = k + 1;
    end

% Allocate room for the persistent spline coefficients
function condition = uniqueSplines( condition )
    [fu, nsfu] = repSplines( condition.sym.fu, 0 );
    
    [dfudp, nsdfudp] = repSplines( condition.sym.dfudp, nsfu );
    condition.sym.dfudp = arMyStr2Sym(dfudp);
    
    [fu, nsfuI] = repInput( fu, nsfu + nsdfudp );

    condition.sym.fu = arMyStr2Sym(fu);    
    condition.splines = zeros( nsfu + nsdfudp + nsfuI, 1 );

% Find the used splines. Replace them with fast splines and store the
% number of spline coefficients required.
% This function makes the following replacement:
%   [mono/pos]splineN(double t, double t1, double p1 ... => [mono/pos]fastsplineN(double t, #, udata_splines__, double t1, double p1
% where udata_splines__ is later replaced by double **splineCache, int *idCache (after symbolic parsing) and # contains the spline index
function [ strOut, nSplines ] = repSplines( fu, offset )
    str = char(fu);
    [~,loc1,loc2] = regexp(str, 'spline(\d*)(\w*)\((\w*),', 'split');
    
    loc1(end+1) = numel(str);
    strOut = str( 1 : loc1(1) );
    for c = 1 : numel( loc2 )
        strOut = [ strOut(1:end-1) 'fast' strOut(end) ];
        strOut = [ strOut str( loc1(c) + 1 : loc2(c) ) ]; %#ok<AGROW>
        strOut = [ strOut num2str(c+offset-1) ', udata_splines__,' ]; %#ok<AGROW>
        strOut = [ strOut str( loc2(c) + 1 : loc1(c+1) ) ]; %#ok<AGROW>
    end
    nSplines = numel(loc2);

% This function replaces instances of inputFunction
function [ strOut, nInput ] = repInput( fu, offset )
    str = char(fu);
    [~,loc1,loc2] = regexp(str, 'inputFunction(\d*)(\w*)\((\w*),', 'split');
    loc1(end+1) = numel(str);
    strOut = str( 1 : loc1(1) );
    for c = 1 : numel( loc2 )
        strOut = strOut(1:end);
        strOut = [ strOut str( loc1(c) + 1 : loc2(c) ) ]; %#ok<AGROW>
        strOut = [ strOut num2str(c+offset-1) ', udata_splines__,' ]; %#ok<AGROW>
        strOut = [ strOut str( loc2(c) + 1 : loc1(c+1) ) ]; %#ok<AGROW>
    end
    nInput = numel(loc2);
    

    
function checksum = addToCheckSum(str, checksum)
algs = {'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if(nargin<2)
    checksum = java.security.MessageDigest.getInstance(algs{2});
end

if iscell(str) 
    str = [str{:}];
end
if(~isempty(str))
    checksum.update(uint8(str(:)));
end


function checkstr = getCheckStr(checksum)
h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

clear checksum

% write ODE headers
function arWriteHFilesCondition(fid, config, condition)

fprintf(fid, '#ifndef _MY_%s\n', condition.fkt);
fprintf(fid, '#define _MY_%s\n\n', condition.fkt);

fprintf(fid, '#include <cvodes/cvodes.h>\n'); 
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <cvodes/cvodes_sparse.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');
fprintf(fid, '#include <cvodes/cvodes_klu.h>\n');
%fprintf(fid, '#include <cvodes/cvodes_superlumt.h>\n');
fprintf(fid, '#include <sundials/sundials_sparse.h>\n');
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

fprintf(fid, ' void fu_%s(void *user_data, double t);\n', condition.fkt);
fprintf(fid, ' void fsu_%s(void *user_data, double t);\n', condition.fkt);
fprintf(fid, ' void fv_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dvdx_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dvdu_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dvdp_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' int fx_%s(realtype t, N_Vector x, N_Vector xdot, void *user_data);\n', condition.fkt);
fprintf(fid, ' void fxdouble_%s(realtype t, N_Vector x, double *xdot_tmp, void *user_data);\n', condition.fkt);
fprintf(fid, ' void fx0_%s(N_Vector x0, void *user_data);\n', condition.fkt);
% fprintf(fid, ' int dfxdx_%s(int N, realtype t, N_Vector x,', condition.fkt); % sundials 2.4.0
fprintf(fid, ' int dfxdx_%s(long int N, realtype t, N_Vector x,', condition.fkt); % sundials 2.5.0
fprintf(fid, 'N_Vector fx, DlsMat J, void *user_data,');
fprintf(fid, 'N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);\n');
fprintf(fid, ' int dfxdx_out_%s(realtype t, N_Vector x, realtype* J, void *user_data);', condition.fkt); % sundials 2.5.0
fprintf(fid, ' int dfxdx_sparse_%s(realtype t, N_Vector x,', condition.fkt); % sundials 2.6.1 with KLU/SuperLU
fprintf(fid, 'N_Vector fx, SlsMat J, void *user_data,'); %DlsMat for Dense solver
fprintf(fid, 'N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);\n');
if(config.useSensiRHS)
    fprintf(fid, ' int fsx_%s(int Ns, realtype t, N_Vector x, N_Vector xdot,', condition.fkt);
    fprintf(fid, 'int ip, N_Vector sx, N_Vector sxdot, void *user_data,');
    fprintf(fid, 'N_Vector tmp1, N_Vector tmp2);\n');
    fprintf(fid, ' int subfsx_%s(int Ns, realtype t, N_Vector x, N_Vector xdot,', condition.fkt);
    fprintf(fid, 'int ip, N_Vector sx, N_Vector sxdot, void *user_data,');
    fprintf(fid, 'N_Vector tmp1, N_Vector tmp2);\n');
end
fprintf(fid, ' void fsx0_%s(int ip, N_Vector sx0, void *user_data);\n', condition.fkt);
fprintf(fid, ' void subfsx0_%s(int ip, N_Vector sx0, void *user_data);\n', condition.fkt);
fprintf(fid, ' void csv_%s(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dfxdp0_%s(realtype t, N_Vector x, double *dfxdp0, void *user_data);\n\n', condition.fkt);
fprintf(fid, ' void dfxdp_%s(realtype t, N_Vector x, double *dfxdp, void *user_data);\n\n', condition.fkt);
fprintf(fid, ' void fz_%s(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x);\n', condition.fkt);
fprintf(fid, ' void fsz_%s(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx);\n\n', condition.fkt);
fprintf(fid, ' void dfzdx_%s(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x);\n', condition.fkt);
fprintf(fid, '#endif /* _MY_%s */\n', condition.fkt);

fprintf(fid,'\n\n\n');


% Write Condition
function arWriteCFilesCondition(fid, matlab_version, config, model, condition, m, c, timedebug)

arFprintf(2, ' -> writing condition m%i c%i, %s...\n', m, c, model.name);

fprintf(fid, '#include "%s.h"\n',  condition.fkt);
fprintf(fid, '#include <cvodes/cvodes.h>\n');    
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <cvodes/cvodes_sparse.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
%fprintf(fid, '#include <cvodes/cvodes_superlumt.h>\n');
fprintf(fid, '#include <sundials/sundials_sparse.h>\n');
fprintf(fid, '#include <cvodes/cvodes_klu.h>\n');
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

% write const vars
fprintf(fid, '%s\n\n', condition.constVars);

% write fu
fprintf(fid, ' void fu_%s(void *user_data, double t)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fu\\n", t);\n'); 
end
if(~isempty(model.us))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    writeCcode(fid, matlab_version, condition, 'fu');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write fsu
fprintf(fid, ' void fsu_%s(void *user_data, double t)\n{\n', condition.fkt);
if(config.useSensis)
    if(isfield(condition.sym, 'dfudp'))
        if(sum(logical(condition.sym.dfudp(:)~=0))>0)
            fprintf(fid, '  UserData data = (UserData) user_data;\n');
            fprintf(fid, '  double *p = data->p;\n');
            
            writeCcode(fid, matlab_version, condition, 'fsu');
        end
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write v
fprintf(fid, ' void fv_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fv\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    writeCcode(fid, matlab_version, condition, 'fv');
end

fprintf(fid, '\n  return;\n}\n\n\n');


% write dvdx
fprintf(fid, ' void dvdx_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t dvdx\\n", t);\n');
end
if(~isempty(model.xs))
    if(config.useSensis || config.useJacobian)
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
        writeCcode(fid, matlab_version, condition, 'dvdx');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dvdu
fprintf(fid, ' void dvdu_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t dvdu\\n", t);\n');
end
if(~isempty(model.us) && ~isempty(model.xs))
    if(config.useSensis)
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
        writeCcode(fid, matlab_version, condition, 'dvdu');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dvdp
fprintf(fid, ' void dvdp_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
	fprintf(fid, '  printf("%%g \\t dvdp\\n", t);\n');
end
if(~isempty(model.xs))
    if(config.useSensis)
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
        if(~isempty(condition.sym.dfvdp))
            writeCcode(fid, matlab_version, condition, 'dvdp');
        end
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write fx
fprintf(fid, ' int fx_%s(realtype t, N_Vector x, N_Vector xdot, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fx\\n", t);\n');
end
fprintf(fid, '  UserData data = (UserData) user_data;\n');
if(~isempty(model.xs))
    fprintf(fid, '  int is;\n');
    fprintf(fid, '  double *qpositivex = data->qpositivex;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *v = data->v;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    fprintf(fid, '  double *xdot_tmp = N_VGetArrayPointer(xdot);\n');
    fprintf(fid, '  fu_%s(data, t);\n', condition.fkt);
    fprintf(fid, '  fv_%s(t, x, data);\n', condition.fkt);
    writeCcode(fid, matlab_version, condition, 'fx');
    fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs));
    fprintf(fid, '    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;\n');
    fprintf(fid, '    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];\n');
    fprintf(fid, '  }\n');
end

fprintf(fid, '\n  return(*(data->abort));\n}\n\n\n');

% write fxdouble
fprintf(fid, ' void fxdouble_%s(realtype t, N_Vector x, double *xdot_tmp, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fxdouble\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  int is;\n');
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *v = data->v;\n');
    fprintf(fid, '  fu_%s(data, t);\n', condition.fkt);
    fprintf(fid, '  fv_%s(t, x, data);\n', condition.fkt);
    writeCcode(fid, matlab_version, condition, 'fx');
    fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs));
    fprintf(fid, '    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;\n');
    fprintf(fid, '  }\n');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write fx0
fprintf(fid, ' void fx0_%s(N_Vector x0, void *user_data)\n{\n', condition.fkt);
if(~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x0_tmp = N_VGetArrayPointer(x0);\n');
    writeCcode(fid, matlab_version, condition, 'fx0');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dfxdx
% fprintf(fid, ' int dfxdx_%s(int N, realtype t, N_Vector x, \n', condition.fkt); % sundials 2.4.0
fprintf(fid, ' int dfxdx_%s(long int N, realtype t, N_Vector x, \n', condition.fkt); % sundials 2.5.0
fprintf(fid, '  \tN_Vector fx, DlsMat J, void *user_data, \n');
fprintf(fid, '  \tN_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\n{\n');
if(timedebug)
    fprintf(fid, '  printf("%%g \\t dfxdx\\n", t);\n');
end

if(~isempty(model.xs))
    if(config.useSensis || config.useJacobian)
        fprintf(fid, '  int is;\n');
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *dvdx = data->dvdx;\n');
        % fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
        fprintf(fid, '  dvdx_%s(t, x, data);\n', condition.fkt);
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs)^2);
        fprintf(fid, '    J->data[is] = 0.0;\n');
        fprintf(fid, '  }\n');
        writeCcode(fid, matlab_version, condition, 'dfxdx');
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs)^2);
        fprintf(fid, '    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;\n');
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return(0);\n}\n\n\n');

% Write dfxdx output function
fprintf(fid, ' int dfxdx_out_%s(realtype t, N_Vector x, realtype* J, void *user_data) {\n', condition.fkt); % sundials 2.5.0
if(~isempty(model.xs))
    if(config.useSensis || config.useJacobian)
        fprintf(fid, '  int is;\n');
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *dvdx = data->dvdx;\n');
        fprintf(fid, '  dvdx_%s(t, x, data);\n', condition.fkt);
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs)^2);
        fprintf(fid, '    J[is] = 0.0;\n');
        fprintf(fid, '  }\n');
        writeCcode(fid, matlab_version, condition, 'dfxdx_out');
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs)^2);
        fprintf(fid, '    if(mxIsNaN(J[is])) J[is] = 0.0;\n');
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return(0);\n}\n\n\n');

% write sparse dfxdx SPARSE (KLU)
fprintf(fid, ' int dfxdx_sparse_%s(realtype t, N_Vector x, \n', condition.fkt); % sundials 2.6.1 with KLU
fprintf(fid, '  \tN_Vector fx, SlsMat J, void *user_data, \n');
fprintf(fid, '  \tN_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\n{\n');
if(timedebug)
    fprintf(fid, '  printf("%%g \\t dfxdx\\n", t);\n');
end

if(~isempty(model.xs))
    if(config.useSensis || config.useJacobian)
        fprintf(fid, '  int is;\n');
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *dvdx = data->dvdx;\n');
        fprintf(fid, '  dvdx_%s(t, x, data);\n', condition.fkt);
        
%         fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(condition.dfxdx_rowVals));
%         fprintf(fid, '    J->data[is] = 0.0;\n');
%         fprintf(fid, '  }\n');         
        fprintf(fid, '  SlsSetToZero(J);\n');
        for j=1:length(condition.dfxdx_rowVals)
            fprintf(fid, '    J->rowvals[%i] = %i', j-1, condition.dfxdx_rowVals(j)); 
            fprintf(fid, ';\n'); 
        end
        fprintf(fid, '\n');
        for j=1:length(condition.dfxdx_colptrs)
            fprintf(fid, '    J->colptrs[%i] = %i', j-1, condition.dfxdx_colptrs(j)); 
            fprintf(fid, ';\n'); 
        end
        fprintf(fid, '\n');
        writeCcode(fid, matlab_version, condition, 'dfxdx_sparse');
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(condition.dfxdx_rowVals));
        fprintf(fid, '    if(mxIsNaN(J->data[is])) J->data[is] = RCONST(0.0);\n');
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return(0);\n}\n\n\n');

% write fsv & fsx
if(config.useSensiRHS)   
    fprintf(fid, ' int fsx_%s(int Ns, realtype t, N_Vector x, N_Vector xdot, \n', condition.fkt);
    fprintf(fid, '  \tint ip, N_Vector sx, N_Vector sxdot, void *user_data, \n');
    fprintf(fid, '  \tN_Vector tmp1, N_Vector tmp2)\n{\n');
    
    if(~isempty(model.xs))
        if(config.useSensis)
            fprintf(fid, '  int is;\n');
            fprintf(fid, '  UserData data = (UserData) user_data;\n');
            fprintf(fid, '  double *p = data->p;\n');
            fprintf(fid, '  double *u = data->u;\n');
            fprintf(fid, '  double *sv = data->sv;\n');
            fprintf(fid, '  double *dvdx = data->dvdx;\n');
            fprintf(fid, '  double *dvdu = data->dvdu;\n');
            fprintf(fid, '  double *dvdp = data->dvdp;\n');
            fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
            
            fprintf(fid, '  double *su = data->su;\n');
            fprintf(fid, '  double *sx_tmp = N_VGetArrayPointer(sx);\n');
            fprintf(fid, '  double *sxdot_tmp = N_VGetArrayPointer(sxdot);\n');
            
            if(timedebug)
                fprintf(fid, '  printf("%%g \\t fsx%%i\\n", t, ip);\n');
            end
            fprintf(fid, '  fsu_%s(data, t);\n', condition.fkt);
            fprintf(fid, '  dvdx_%s(t, x, data);\n', condition.fkt);
            fprintf(fid, '  dvdu_%s(t, x, data);\n', condition.fkt);
            fprintf(fid, '  dvdp_%s(t, x, data);\n', condition.fkt);
            
            % Flux sensitivities
            fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(condition.sv));
            fprintf(fid, '    sv[is] = 0.0;\n');
            fprintf(fid, '  }\n');
            
            writeCcode(fid, matlab_version, condition, 'fsv1');
            fprintf(fid, '  switch (ip) {\n');
            for j2=1:size(condition.sym.dvdp,2)
                fprintf(fid, '    case %i: {\n', j2-1);
                writeCcode(fid, matlab_version, condition, 'fsv2', j2);
                fprintf(fid, '    } break;\n');
            end
            fprintf(fid, '  }\n');
            writeCcode(fid, matlab_version, condition, 'fsx');
            
            % Add sensitivity RHS contributions corresponding to the compartment volumes
            if ( isfield( condition.sym, 'dfcdp' ) )
                for svs = 1 : length(condition.sym.fsx)
                    sxdot_tmp{svs} = sprintf( 'sxdot_tmp[%d]', svs ); %#ok
                end
                sxdot_tmp = arMyStr2Sym(repmat(sxdot_tmp, size(condition.sym.dvdp,2), 1));
                condition.sym.dfcdp2 = (sxdot_tmp.' + condition.sym.dfcdp);
                
                fprintf(fid, '  switch (ip) {\n');
                for j2=1:size(condition.sym.dvdp,2)
                    fprintf(fid, '    case %i: {\n', j2-1);
                    writeCcode(fid, matlab_version, condition, 'dfcdp2', j2);
                    fprintf(fid, '    } break;\n');
                end
                fprintf(fid, '  }\n');
            end
            
            fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs));
            fprintf(fid, '    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;\n');
            fprintf(fid, '  }\n');
        end
    end
    fprintf(fid, '\n  return(0);\n}\n\n\n');
    
    % Function which can be called to get a smaller sensitivity matrix if
    % not all sensitivities are required.
    fprintf(fid, ' int subfsx_%s(int Ns, realtype t, N_Vector x, N_Vector xdot, \n', condition.fkt);
    fprintf(fid, '  \tint ip, N_Vector sx, N_Vector sxdot, void *user_data, \n');
    fprintf(fid, '  \tN_Vector tmp1, N_Vector tmp2)\n {\n');
    fprintf(fid, '   UserData data = (UserData) user_data;\n');   
    fprintf(fid, '   return fsx_%s(Ns, t, x, xdot, data->sensIndices[ip], sx, sxdot, user_data, tmp1, tmp2);\n', condition.fkt);
    %fprintf(fid, ' return 0;');
    fprintf(fid, ' };\n\n');    
end

% compute flux sensitivity (only for output)
fprintf(fid, ' void csv_%s(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data)\n{\n', condition.fkt);
if(~isempty(model.xs))
    if(config.useSensis)
        fprintf(fid, '  int is;\n');
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *sv = data->sv;\n');
        fprintf(fid, '  double *dvdx = data->dvdx;\n');
        fprintf(fid, '  double *dvdu = data->dvdu;\n');
        fprintf(fid, '  double *dvdp = data->dvdp;\n');
        fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');

        fprintf(fid, '  double *su = data->su;\n');
        fprintf(fid, '  double *sx_tmp = N_VGetArrayPointer(sx);\n');

        % Flux sensitivities
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(condition.sv));
        fprintf(fid, '    sv[is] = 0.0;\n');
        fprintf(fid, '  }\n');

        writeCcode(fid, matlab_version, condition, 'fsv1');
        fprintf(fid, '  switch (ip) {\n');
        for j2=1:size(condition.sym.dvdp,2)
            fprintf(fid, '    case %i: {\n', j2-1);
            writeCcode(fid, matlab_version, condition, 'fsv2', j2);
            fprintf(fid, '    } break;\n');
        end
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');    

% write fsx0
if (ispc)
    % On Windows SDK 7.1; code optimization sometimes leads to a rare segfault
    % using specific expressions. I am currently investigating this, but
    % this is to prevent it from breaking my current models.
    fprintf(fid, ' #pragma optimize("", off)\n');
end

fprintf(fid, ' void fsx0_%s(int ip, N_Vector sx0, void *user_data)\n{\n', condition.fkt);
if(~isempty(model.xs))
    if(config.useSensis)
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *sx0_tmp = N_VGetArrayPointer(sx0);\n');
        
        % Equations
        fprintf(fid, '  switch (ip) {\n');
        for j2=1:size(condition.sym.fsx0,2)
            if(sum(logical(condition.sym.fsx0(:,j2)~=0))>0)
                fprintf(fid, '    case %i: {\n', j2-1);
                writeCcode(fid, matlab_version, condition, 'fsx0', j2);
                fprintf(fid, '    } break;\n');
            end
        end
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% Function which can be called to get a smaller sensitivity matrix if
% not all sensitivities are required.
fprintf(fid, ' void subfsx0_%s(int ip, N_Vector sx0, void *user_data)\n {\n', condition.fkt);
fprintf(fid, '   UserData data = (UserData) user_data;\n');
fprintf(fid, '   fsx0_%s(data->sensIndices[ip], sx0, user_data);\n', condition.fkt);
fprintf(fid, ' };\n\n');

if (ispc)
    fprintf(fid, ' #pragma optimize("", on)\n');
end

% write dfxdp0 (total derivative at initial condition)
fprintf(fid, ' void dfxdp0_%s(realtype t, N_Vector x, double *dfxdp0, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
	fprintf(fid, '  printf("%%g \\t dfxdp0\\n", t);\n');
end
if(~isempty(model.xs))
    if(config.useSensis)
        if(~isempty(condition.sym.dfxdp0))
            fprintf(fid, '  int is;\n');
            fprintf(fid, '  UserData data = (UserData) user_data;\n');
            fprintf(fid, '  double *p = data->p;\n');
            fprintf(fid, '  double *u = data->u;\n');
            fprintf(fid, '  double *dvdp = data->dvdp;\n');
            fprintf(fid, '  double *dvdx = data->dvdx;\n');
            fprintf(fid, '  double *dvdu = data->dvdu;\n');
            fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
            writeCcode(fid, matlab_version, condition, 'dfxdp0');
            fprintf(fid, '  for (is=0; is<%i; is++) {\n', numel(condition.sym.dfxdp0));
            fprintf(fid, '    if(mxIsNaN(dfxdp0[is])) dfxdp0[is] = 0.0;\n');
            fprintf(fid, '  }\n');
        end
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dfxdp (partial derivative w.r.t. p and keeping x constant)
fprintf(fid, ' void dfxdp_%s(realtype t, N_Vector x, double *dfxdp, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
	fprintf(fid, '  printf("%%g \\t dfxdp\\n", t);\n');
end
if(~isempty(model.xs))
    if(config.useSensis)
        if(~isempty(condition.sym.dfxdp))
            fprintf(fid, '  int is;\n');
            fprintf(fid, '  UserData data = (UserData) user_data;\n');
            fprintf(fid, '  double *p = data->p;\n');
            fprintf(fid, '  double *u = data->u;\n');
            fprintf(fid, '  double *dvdp = data->dvdp;\n');
            fprintf(fid, '  double *dvdx = data->dvdx;\n');
            fprintf(fid, '  double *dvdu = data->dvdu;\n');
            fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
            writeCcode(fid, matlab_version, condition, 'dfxdp');
            fprintf(fid, '  for (is=0; is<%i; is++) {\n', numel(condition.sym.dfxdp));
            fprintf(fid, '    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;\n');
            fprintf(fid, '  }\n');
        end
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write z
fprintf(fid, ' void fz_%s(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x){\n', condition.fkt);
if(~isempty(model.zs))
    writeCcode(fid, matlab_version, condition, 'fz');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write sz
fprintf(fid, ' void fsz_%s(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx){\n', condition.fkt);
if(~isempty(model.zs))
    if(config.useSensis)
        fprintf(fid, '  int jp;\n');
        fprintf(fid, '  for (jp=0; jp<np; jp++) {\n');
        writeCcode(fid, matlab_version, condition, 'fsz1');
        fprintf(fid, '  };\n');
        fprintf(fid, '\n');
        writeCcode(fid, matlab_version, condition, 'fsz2');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dfzdx
fprintf(fid, ' void dfzdx_%s(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdxs, double *z, double *p, double *u, double *x){\n', condition.fkt);
if(~isempty(model.zs))
    if(config.useSensis)
        writeCcode(fid, matlab_version, condition, 'dfzdx');        
        fprintf(fid, '\n');        
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');


% write data headers
function arWriteHFilesData(fid, data)

fprintf(fid, '#ifndef _MY_%s\n', data.fkt);
fprintf(fid, '#define _MY_%s\n\n', data.fkt);

fprintf(fid, '#include <cvodes/cvodes.h>\n'); 
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <cvodes/cvodes_sparse.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
%fprintf(fid, '#include <cvodes/cvodes_superlumt.h>\n');
fprintf(fid, '#include <sundials/sundials_sparse.h>\n');
fprintf(fid, '#include <cvodes/cvodes_klu.h>\n');
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

fprintf(fid, ' void fy_%s(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z);\n', data.fkt);
fprintf(fid, ' void fystd_%s(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z);\n', data.fkt);
fprintf(fid, ' void fsy_%s(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz);\n', data.fkt);
fprintf(fid, ' void fsystd_%s(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz);\n\n', data.fkt);
fprintf(fid, ' void fy_scale_%s(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx);\n', data.fkt);
fprintf(fid, '#endif /* _MY_%s */\n', data.fkt);
fprintf(fid,'\n\n\n');


% Write Data
function arWriteCFilesData(fid, matlab_version, config, m, c, d, data)

arFprintf(2, ' -> writing data m%i d%i -> c%i, %s...\n', m, d, c, data.name);

fprintf(fid, '#include "%s.h"\n',  data.fkt);
fprintf(fid, '#include <cvodes/cvodes.h>\n');    
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <cvodes/cvodes_sparse.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
%fprintf(fid, '#include <cvodes/cvodes_superlumt.h>\n');
fprintf(fid, '#include <sundials/sundials_sparse.h>\n');
fprintf(fid, '#include <cvodes/cvodes_klu.h>\n');
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

% write const vars
fprintf(fid, '%s\n\n', data.constVars);

% write y
fprintf(fid, ' void fy_%s(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z){\n', data.fkt);
writeCcode(fid, matlab_version, data, 'fy');
fprintf(fid, '\n  return;\n}\n\n\n');

% write ystd
fprintf(fid, ' void fystd_%s(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z){\n', data.fkt);
writeCcode(fid, matlab_version, data, 'fystd');
fprintf(fid, '\n  return;\n}\n\n\n');

% write sy
fprintf(fid, ' void fsy_%s(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz){\n', data.fkt);
if(config.useSensis)
    writeCcode(fid, matlab_version, data, 'fsy');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write systd
fprintf(fid, ' void fsystd_%s(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz){\n', data.fkt);
if(config.useSensis)
    writeCcode(fid, matlab_version, data, 'fsystd');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write y_scale
fprintf(fid, ' void fy_scale_%s(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx){\n', data.fkt);
if(~isempty(data.sym.y_scale))
	writeCcode(fid, matlab_version, data, 'y_scale');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write C code
function writeCcode(fid, matlab_version, cond_data, svar, ip)
    
if(strcmp(svar,'fv'))
    cstr = ccode2(cond_data.sym.fv(:), matlab_version);
    cvar =  'data->v';
elseif(strcmp(svar,'dvdx'))
    cstr = ccode2(cond_data.sym.dfvdx(:), matlab_version);
    cvar =  'data->dvdx';
elseif(strcmp(svar,'dvdu'))
    cstr = ccode2(cond_data.sym.dfvdu(:), matlab_version);
    cvar =  'data->dvdu';
elseif(strcmp(svar,'dvdp'))
    cstr = ccode2(cond_data.sym.dfvdp(:), matlab_version);
    cvar =  'data->dvdp';
elseif(strcmp(svar,'fx'))
    cstr = ccode2(cond_data.sym.fx(:), matlab_version);
    for j=find(cond_data.sym.fx(:)' == 0)
        cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
    end
    cvar =  'xdot_tmp';
elseif(strcmp(svar,'fx0'))
    cstr = ccode2(cond_data.sym.fpx0(:), matlab_version);
    cvar =  'x0_tmp';
elseif(strcmp(svar,'dfxdx'))
    cstr = ccode2(cond_data.sym.dfxdx(:), matlab_version);
%     for j=find(cond_data.sym.dfxdx(:)' == 0)
%         cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
%     end
    cvar =  'J->data';
elseif(strcmp(svar,'dfxdx_sparse'))
    cstr = ccode2(cond_data.sym.dfxdx_nonzero(:), matlab_version);    
    cvar =  'J->data';
elseif(strcmp(svar,'dfxdx_out'))
    cstr = ccode2(cond_data.sym.dfxdx(:), matlab_version);
    cvar =  'J';
elseif(strcmp(svar,'fsv1'))
    cstr = ccode2(cond_data.sym.fsv1, matlab_version);
    cvar =  'sv';
elseif(strcmp(svar,'fsv2'))
    cstr = ccode2(cond_data.sym.dvdp(:,ip), matlab_version);
    cvar =  '    sv';
elseif(strcmp(svar,'fsx'))
    cstr = ccode2(cond_data.sym.fsx, matlab_version);
    for j=find(cond_data.sym.fsx' == 0)
        cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
    end
    cvar =  'sxdot_tmp';
elseif(strcmp(svar,'dfcdp2'))
    cstr = ccode2(cond_data.sym.dfcdp2(:,ip), matlab_version);
    cvar =  'sxdot_tmp';    
elseif(strcmp(svar,'fsx0'))
    cstr = ccode2(cond_data.sym.fsx0(:,ip), matlab_version);
    cvar =  '    sx0_tmp';
elseif(strcmp(svar,'fu'))
    cstr = ccode2(cond_data.sym.fu(:), matlab_version);
    cvar =  'data->u';
elseif(strcmp(svar,'fsu'))
    cstr = ccode2(cond_data.sym.dfudp(:), matlab_version);
    cvar =  'data->su';
elseif(strcmp(svar,'fz'))
    cstr = ccode2(cond_data.sym.fz(:), matlab_version);
    cvar =  'z';
elseif(strcmp(svar,'dfzdx'))
    cstr = ccode2(cond_data.sym.dfzdx(:), matlab_version);
    cvar =  '    dfzdxs';
elseif(strcmp(svar,'fsz1'))
    cstr = ccode2(cond_data.sym.fsz1, matlab_version);
    for j=find(cond_data.sym.fsz1' == 0)
        cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
    end
    cvar =  '    sz';
elseif(strcmp(svar,'fsz2'))
    cstr = ccode2(cond_data.sym.fsz2(:), matlab_version);
    cvar =  'sz';
elseif(strcmp(svar,'fy'))
    cstr = ccode2(cond_data.sym.fy(:), matlab_version);
    cvar =  'y';
elseif(strcmp(svar,'y_scale'))
    cstr = ccode2(cond_data.sym.y_scale(:), matlab_version);
    cvar =  'y_scale';
elseif(strcmp(svar,'fystd'))
    cstr = ccode2(cond_data.sym.fystd(:), matlab_version);
    cvar =  'ystd';
elseif(strcmp(svar,'fsy'))
    cstr = ccode2(cond_data.sym.fsy(:), matlab_version);
    cvar =  'sy';
elseif(strcmp(svar,'fsystd'))
    cstr = ccode2(cond_data.sym.fsystd(:), matlab_version);
    cvar =  'systd';
elseif(strcmp(svar,'dfxdp0'))
    cstr = ccode2(cond_data.sym.dfxdp0(:), matlab_version);
    cvar =  'dfxdp0';
elseif(strcmp(svar,'dfxdp'))
    cstr = ccode2(cond_data.sym.dfxdp(:), matlab_version);
    cvar =  'dfxdp';
else
    error('unknown %s', svar);
end

cstr = strrep(cstr, 't0', [cvar '[0]']);
cstr = strrep(cstr, '][0]', ']');
cstr = strrep(cstr, 'T[', [cvar '[']);
cstr = strrep(cstr, 'udata_splines__', 'data->splines, data->splineIndices');

% % debug
% fprintf('\n\n');
% if(config.isMaple)
%     for j=1:length(cstr)
%         fprintf('%s\n', cstr{j});
%     end
% else
%     fprintf('%s', cstr);
% end
% fprintf('\n');

% Find instances of heaviside functions (these need special treatment)
%if (exist('ar','var') && ~isempty(ar) && isfield( ar.config, 'accurateSteps' ) )
%    if ( ar.config.accurateSteps == 1 )
%        cstr = replaceWithinFunc(cstr, 'heaviside', 't', 'data->t');
%    end
%end

if(~(length(cstr)==1 && isempty(cstr{1})))
    if(strcmp(svar,'fy'))
        cstr = strrep(cstr, 'x[', 'x[nx*ntlink*iruns+itlink+ntlink*');
        cstr = strrep(cstr, 'z[', 'z[nz*ntlink*iruns+itlink+ntlink*');
        cstr = strrep(cstr, 'y[', 'y[ny*nt*iruns+it+nt*');
        cstr = strrep(cstr, 'u[', 'u[itlink+ntlink*');
    elseif(strcmp(svar,'y_scale'))
        cstr = strrep(cstr, 'x[', 'x[nx*ntlink*iruns+itlink+ntlink*');
        cstr = strrep(cstr, 'z[', 'z[nz*ntlink*iruns+itlink+ntlink*');
        cstr = strrep(cstr, 'y_scale[', 'y_scale[ny*nt*iruns+it+nt*');
        cstr = strrep(cstr, 'u[', 'u[itlink+ntlink*');
    elseif(strcmp(svar,'fsy'))
        cstr = strrep(cstr, 'x[', 'x[itlink+ntlink*');
        cstr = strrep(cstr, 'z[', 'z[itlink+ntlink*');
        cstr = strrep(cstr, 'u[', 'u[itlink+ntlink*');
        cstr = strrep(cstr, 'y[', 'y[it+nt*');        
    elseif(strcmp(svar,'fystd'))
        cstr = strrep(cstr, 'x[', 'x[itlink+ntlink*');
        cstr = strrep(cstr, 'z[', 'z[itlink+ntlink*');
        cstr = strrep(cstr, 'u[', 'u[itlink+ntlink*');
        cstr = strrep(cstr, 'y[', 'y[it+nt*');
        cstr = strrep(cstr, 'ystd[', 'ystd[it+nt*');
    elseif(strcmp(svar,'fsystd'))
        cstr = strrep(cstr, 'x[', 'x[itlink+ntlink*');
        cstr = strrep(cstr, 'z[', 'z[itlink+ntlink*');
        cstr = strrep(cstr, 'u[', 'u[itlink+ntlink*');
        cstr = strrep(cstr, 'y[', 'y[it+nt*');
        cstr = strrep(cstr, 'ystd[', 'ystd[it+nt*');
        
    elseif(strcmp(svar,'fz'))
        cstr = strrep(cstr, 'x[', 'x[nx*nt*iruns+it+nt*');
        cstr = strrep(cstr, 'z[', 'z[nz*nt*iruns+it+nt*');
        cstr = strrep(cstr, 'u[', 'u[nu*nt*iruns+it+nt*');
    elseif(strcmp(svar,'dfzdx'))
        cstr = strrep(cstr, 'x[', 'x[nx*nt*iruns+it+nt*');
        cstr = strrep(cstr, 'dfzdxs[', 'dfzdxs[nx*nt*iruns+it+nt*');
        cstr = strrep(cstr, 'u[', 'u[nu*nt*iruns+it+nt*');
    elseif(strcmp(svar,'fsz1'))
        cstr = strrep(cstr, 'x[', 'x[it+nt*');
        cstr = strrep(cstr, 'z[', 'z[it+nt*');
        cstr = strrep(cstr, 'u[', 'u[it+nt*');
        cstr = strrep(cstr, 'sx[it+nt*', sprintf('sx[it + nt*%i*jp + nt*', length(cond_data.sx)));
        cstr = strrep(cstr, 'sz[it+nt*', sprintf('sz[it + nt*%i*jp + nt*', length(cond_data.sz)));
        cstr = strrep(cstr, 'su[it+nt*', sprintf('su[it + nt*%i*jp + nt*', length(cond_data.su)));
    elseif(strcmp(svar,'fsz2'))
        cstr = strrep(cstr, 'x[', 'x[it+nt*');
        cstr = strrep(cstr, 'z[', 'z[it+nt*');
        cstr = strrep(cstr, 'u[', 'u[it+nt*');
        cstr = strrep(cstr, '=', '+=');
        
    elseif(strcmp(svar,'fsv1'))
        cstr = strrep(cstr, 'su[', sprintf('su[(ip*%i)+',length(cond_data.su)));
        cstr = strrep(cstr, 'x[', 'x_tmp[');
		cstr = strrep(cstr, 'dvdx_tmp', 'dvdx');
    elseif(strcmp(svar,'fsv2'))
        cstr = strrep(cstr, '=', '+=');
        cstr = strrep(cstr, 'x[', 'x_tmp[');
		cstr = strrep(cstr, 'dvdx_tmp', 'dvdx');
        
    else
        cstr = strrep(cstr, 'x[', 'x_tmp[');
		cstr = strrep(cstr, 'dvdx_tmp', 'dvdx');
    end
end

fprintf(fid, '%s\n', cstr);

% % debug
% fprintf('\n\n%s\n', cstr);


function writeSimuCalcFunctions(debug_mode)

global ar

% Functions
fid = fopen(['./Compiled/' ar.info.c_version_code '/arSimuCalcFunctions.c'], 'W');

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(~debug_mode)
            fprintf(fid, '#include "%s.h"\n', ar.model(m).condition(c).fkt);
        else
            fprintf(fid, '#include "%s.c"\n', ar.model(m).condition(c).fkt);
        end
    end
    
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            if(~debug_mode)
                fprintf(fid, '#include "%s.h"\n', ar.model(m).data(d).fkt);
            else
                fprintf(fid, '#include "%s.c"\n', ar.model(m).data(d).fkt);
            end
        end
    end
end
fprintf(fid, '\n');

% map CVodeInit to fx
fprintf(fid, ' int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) return CVodeInit(cvode_mem, fx_%s, RCONST(t), x);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '  return(-1);\n');
fprintf(fid, '}\n\n');

% map fx
fprintf(fid, ' void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) fxdouble_%s(t, x, xdot, user_data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map fx0
fprintf(fid, ' void fx0(N_Vector x0, void *user_data, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) fx0_%s(x0, data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map CVDlsSetDenseJacFn to dfxdx
fprintf(fid, ' int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic, int setSparse){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i) & (setSparse==0)){ \n return CVDlsSetDenseJacFn(cvode_mem, dfxdx_%s);\n', ...   
            m-1, c-1, ar.model(m).condition(c).fkt);
        fprintf(fid, ' \n }else if((im==%i) & (ic==%i) & (setSparse==1)){ \n return CVSlsSetSparseJacFn(cvode_mem, dfxdx_sparse_%s);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
        fprintf(fid, '\n}\n');
    end
end
fprintf(fid, '  return(-1);\n');
fprintf(fid, '}\n\n');

% map dfxdx output function
fprintf(fid, ' void getdfxdx(int im, int ic, realtype t, N_Vector x, realtype *J, void *user_data){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) { dfxdx_out_%s(t, x, J, user_data); return; }\n', ...   
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '\n}\n\n');

% map fsx0
fprintf(fid, ' void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic, int sensitivitySubset) {\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
fprintf(fid, '  if ( sensitivitySubset == 0 ) {\n' );
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '    if((im==%i) & (ic==%i)) fsx0_%s(is, sx_is, data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '  } else {\n' );
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '    if((im==%i) & (ic==%i)) subfsx0_%s(is, sx_is, data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '  }\n' );
fprintf(fid, '}\n\n');

% map flux sensitivities <== csv is here
fprintf(fid, ' void csv(realtype t, N_Vector x, int ip, N_Vector sx, void *user_data, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) csv_%s(t, x, ip, sx, data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map CVodeSensInit1 to fsx
fprintf(fid, ' int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic, int sensitivitySubset){\n');
if(ar.config.useSensiRHS)
    fprintf(fid, '  if (sensirhs == 1) {\n');
    fprintf(fid, '    if (sensitivitySubset == 0) {\n');
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            fprintf(fid, '      if((im==%i) & (ic==%i)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_%s, sx);\n', ...
                m-1, c-1, ar.model(m).condition(c).fkt);
        end
    end
    fprintf(fid, '    } else {\n');
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            fprintf(fid, '      if((im==%i) & (ic==%i)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, subfsx_%s, sx);\n', ...
                m-1, c-1, ar.model(m).condition(c).fkt);
        end
    end
    fprintf(fid, '  }\n');
    fprintf(fid, '  } else {\n');
end
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '    if((im==%i) & (ic==%i)) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);\n', ...
            m-1, c-1);
    end
end
if(ar.config.useSensiRHS)
    fprintf(fid, '  }\n');
end
fprintf(fid, '  return(-1);\n');
fprintf(fid, '}\n\n');

% map fu
fprintf(fid, ' void fu(void *user_data, double t, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) fu_%s(data, t);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fsu
fprintf(fid, ' void fsu(void *user_data, double t, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) fsu_%s(data, t);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fv
fprintf(fid, ' void fv(void *user_data, double t, N_Vector x, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) fv_%s(t, x, data);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fsv
fprintf(fid, ' void fsv(void *user_data, double t, N_Vector x, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) {\n\tdvdp_%s(t, x, data);\n\tdvdu_%s(t, x, data);\n\tdvdx_%s(t, x, data);\n}\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt, ar.model(m).condition(c).fkt, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map dfxdp0
fprintf(fid, ' void dfxdp0(void *user_data, double t, N_Vector x, double *dfxdp0, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) dfxdp0_%s(t, x, dfxdp0, data);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map dfxdp
fprintf(fid, ' void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) dfxdp_%s(t, x, dfxdp, data);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fz
fprintf(fid, 'void fz(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *z, double *p, double *u, double *x, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) fz_%s(t, nt, it, nz, nx, nu, iruns, z, p, u, x);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map dfzdx
fprintf(fid, 'void dfzdx(double t, int nt, int it, int nz, int nx, int nu, int iruns, double *dfzdx, double *z, double *p, double *u, double *x, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) dfzdx_%s(t, nt, it, nz, nx, nu, iruns, dfzdx, z, p, u, x);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map fsz
fprintf(fid, 'void fsz(double t, int nt, int it, int np, double *sz, double *p, double *u, double *x, double *z, double *su, double *sx, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if((im==%i) & (ic==%i)) fsz_%s(t, nt, it, np, sz, p, u, x, z, su, sx);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map fy
fprintf(fid, ' void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y, double *p, double *u, double *x, double *z, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if((im==%i) & (id==%i)) fy_%s(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y, p, u, x, z);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fy_scale
fprintf(fid, ' void fy_scale(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int nz, int iruns, double *y_scale, double *p, double *u, double *x, double *z, double *dfzdx, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if((im==%i) & (id==%i)) fy_scale_%s(t, nt, it, ntlink, itlink, ny, nx, nz, iruns, y_scale, p, u, x, z, dfzdx);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fystd
fprintf(fid, ' void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, double *z, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if((im==%i) & (id==%i)) fystd_%s(t, nt, it, ntlink, itlink, ystd, y, p, u, x, z);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fsy
fprintf(fid, ' void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *z, double *su, double *sx, double *sz, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if((im==%i) & (id==%i)) fsy_%s(t, nt, it, ntlink, itlink, sy, p, u, x, z, su, sx, sz);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fsystd
fprintf(fid, ' void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *z, double *sy, double *su, double *sx, double *sz, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if((im==%i) & (id==%i)) fsystd_%s(t, nt, it, ntlink, itlink, systd, p, y, u, x, z, sy, su, sx, sz);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% for arSSACalc
% call to fu and fv
fprintf(fid, '/* for arSSACalc.c */\n\n');
fprintf(fid, ' void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if((im==%i) & (ic==%i)) {\n', m-1, c-1);
        fprintf(fid, '    fu_%s(data, t);\n', ar.model(m).condition(c).fkt);
        fprintf(fid, '    fv_%s(t, x, data);\n', ar.model(m).condition(c).fkt);
        fprintf(fid, '  }\n');
	end
end
fprintf(fid, '}\n\n');

fclose(fid);


function J = myJacobian(F,x)
% function checks if first argument is empty to provide R2013b compatibility.
% If F is not empty, the built-in jacobian is called. Else, the function returns 
% an empty sym in the right dimensions.

if(~isempty(F))
    if(~isempty(x))
        J = jacobian(F,x);
    else
        J = arMyStr2Sym(NaN(length(F),0));
    end
else
    J = arMyStr2Sym(NaN(0,length(x)));
end

% Replace special functions before converting to symbolic expression
function s = mySym( s, specialFunc )
    if ( isempty( specialFunc ) )
        s = arMyStr2Sym(s);
    else
        for a = 1 : size( s, 1 )
            for b = 1 : size( s, 2 )
                s{a,b} = replaceFunctions( s{a,b}, specialFunc, 1 );
            end
        end
        s = arMyStr2Sym(s);
    end
    
% convert sym array to string array
function a = sym2str(b)
    nonzeros = logical(b~=0);
    a = cell(size(b));
    for j=1:length(b)
        if ( nonzeros(j) )
            a{j} = char(b(j));
        else
            a{j} = '0';
        end
    end
    
% convert sym array to string array
function a = sym2str_old(b)
    a = cell(size(b));
    for j=1:length(b)
        a{j} = char(b(j));
    end
    
function strG = quickScan( str )
    c = 1;
    depth = 0;
    while ( c < numel( str ) )
        if ( str(c) == '(' )
            depth = depth + 1;
        end
        if ( str(c) == ')' )
            depth = depth - 1;
            if ( depth < 0 )
                strG = str(1:c);
                return
            end
        end
        c = c + 1;
    end
    strG = str;
    
% This function maps matlab symbolic toolbox derivatives to derivatives we
% can use in C
%
%  D([#], name)(args) => Dname(args, %d)
%  diff(name(args), args(#)) => Dname(args, #)
function str = replaceDerivative( str )
    % Pattern that matches the derivatives D([#], func)(args)
    pattern = 'D[\(][\[](\d+)[\]][\,](\s*)(\w*)[\)][\(]';
    
    % Transform D([#], name)(args) => Dname(args, %d)
    % We need to explicitly match up the brackets, since regexps can't do
    % this reliably, we have to loop over it.
    [tokens,startloc,endloc] = regexp( str, pattern, 'tokens' );
    
    % Generate replacement blocks
    oldString = cell( 1, numel( startloc ) );
    to = cell( 1, numel( startloc ) );
    for i = 1 : numel( startloc )
        call = str(startloc(i):endloc(i));
        block = quickScan( str(endloc(i)+1:end) );
        oldString{i} = [ call block ];
        to{i} = sprintf( 'D%s(%s, %s)', tokens{i}{3}, block(1:end-1), tokens{i}{1} );
    end
    
    % Perform the replacements
    for i = 1 : numel( oldString )
        str = strrep( str, oldString{i}, to{i} );
    end
        
    % Pattern which matches the other derivative structure
    pattern2 = 'diff[\(]([\w\(\)]+)[\(]([\[\]\(\)\^\/\*\+\-\.\s,\w\d]*)[\)][,](\s)([\[\]\w]*)[\)]';
        
    % Performs regexprep which transforms diff(name(args), args(#)) => Dname(args, #)
    [~,~,~,total,matches]=regexp(str, pattern2);

    for jm = 1 : numel(matches)
        indep = matches{jm}{4};
        variables = strtrim(strsplit(matches{jm}{2}, ','));

        % Find which variable we're deriving w.r.t. to
        ID = find(strcmp(variables, indep));

        % Write the new string
        fNew = sprintf( 'D%s(%s, %d)', matches{jm}{1}, matches{jm}{2}, ID );

        str = strrep( str, total{jm}, fNew );
    end

% Safely map derivatives to the appropriate C functions
%   pattern replaces D([#], func)(args) to Dfunc(args, floor(#/2)) 
function str = repSplineDer( str )

    % Pattern that matches the derivatives D([#], func)(args)
    pattern = 'D[\(][\[](\d+)[\]][\,]\s([\w\(\)]*)[\)][\(]([\[\]\^\/\*\+\-\.\s,\w\d]*)[\)]';
    
    % Compute the mask for the printf
    % Performs regexprep which transforms D([#], name)(args) => Dname(args, %d)
    mask = regexprep(str, pattern,'D$2($3,%d)');
    
    % Now the derivative IDs are computed with another
    % regexprep (divided by two and floored)
    chunks = regexp(str, pattern, 'match');
    values = regexprep(chunks, pattern, '$1');
    values = floor( cellfun(@str2num, values) / 2 );

    str = sprintf( mask, values );  

function prepareBecauseOfRepeatedCompilation
    global ar
    for m=1:length(ar.model)
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).p = ar.model(m).data(d).pold;
        end
    end
        

function cstr = ccode2(T, matlab_version)

    % If this matrix or value is empty, do not attempt to generate C-code
    if (numel(T) == 0)
        cstr = '';
        return;
    end
    
    try
        T = arMyStr2Sym( replaceDerivative( char(T) ) );
    catch
        char(T)
        error('Failure attempting to replace derivatives in');
    end

    % R2015b compatibility fix
    if(matlab_version>=8.6)
        sym_str = sym2str(T);
        if all(strcmp('0',sym_str))
            cstr = char;
        else
            cstr = ccode(T);
        end
    else
        cstr = ccode(T);
    end
    
    % Fixes bug with ccode for single line syms
    if ( length( cstr ) > 8 )
        if ( strcmp(cstr(1:9), '  _assign') )
            cstr = strrep( [cstr(1:end-2) ';'], '_assign(t0,', 't0=' );
        end     
    end
    
function message = invalidSym( message, symvariable, names )
    if ( exist( 'names', 'var' ) )
        for a = 1 : numel( symvariable )
            try
                arMyStr2Sym(symvariable(a));
            catch
                message = sprintf( '%s\n%s -> %s\n', message, names{a}, symvariable{a} );
            end
        end
    else
        for a = 1 : numel( symvariable )
            try
                arMyStr2Sym(symvariable(a));
            catch
                message = sprintf( '%s\n%s\n', message, symvariable{a} );
            end
        end
    end
    
    
function [ dvdx, dvdu ] = cacheShortTerms( model )
    dvdx = sym(zeros(length(model.vs), length(model.xs)));
    for j=1:length(model.vs)
        for i=1:length(model.xs)
            dvdx(j,i) = arMyStr2Sym(sprintf('dvdx[%i]', j + (i-1)*length(model.vs)));
        end
    end

    dvdu = arMyStr2Sym(zeros(length(model.vs), length(model.us)));
    for j=1:length(model.vs)
        for i=1:length(model.us)
            dvdu(j,i) = arMyStr2Sym(sprintf('dvdu[%i]', j + (i-1)*length(model.vs)));
        end
    end
