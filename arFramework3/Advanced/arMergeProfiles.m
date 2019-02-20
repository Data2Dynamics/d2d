% PLE = arMergeProfiles(directory, [forceLoad], [filterProfile], [appendProfiles])
%
%  merge profiles from various runs
%
%  directory        Should point to results directory. Typically, the profiles 
%                   will be subdirectories in this directory; each with a 
%                   subdirectory /PLE/
%                   Alternatively, an ar struct with profiles can also be
%                   passed to the directory argument.
%  forceLoad        Force reload of the model [false]
%  filterProfile    Smoothen profile? [false]
%
%  PLE              Yields the merged PLE struct as the optional output. If
%                   no output is demanded, merged ple struct is saved to ar.ple
%
%  appendProfiles   Append profiles of the same parameter in different workspaces instead of
%                   replacing them  [false]
%
%  If no model is loaded, an incompatible model is loaded or forceLoad is
%  set to one, arMergeProfiles will attempt to load a model. For this, it
%  will look one directory up from Directory for a mat file named
%  arStruct.mat

function vargout = arMergeProfiles( directory, forceLoad, filterProfile, appendProfiles)

    global ar;
    
    if ( ~exist( 'filterProfile', 'var' ) || isempty(filterProfile))
        filterProfile = 0;
    end
    if ( ~exist( 'appendProfiles', 'var' ) )
        appendProfiles = false;
    end 
    
    if ischar(directory)
        s = dir(directory);

        if ( ~exist( 'forceLoad', 'var' ) || isempty(forceLoad))
            forceLoad = 0;
        end
        if isempty( ar ) || forceLoad
            load( [ directory '/../arStruct.mat' ] );
        end
        
        % Assemble ple struct
        PLEs = ar.ple;
        PLEs.donePars = {};
        for a = 1 : length( s )
            if ( ~strcmp( s(a).name, '.' ) && ~strcmp( s(a).name, '..' ) )
                path = [ directory '/' s(a).name ];
                if ( isdir( path ) )
                    path
                    PLEs = loadSingle( PLEs, path, filterProfile, appendProfiles);
                end
            end
        end
    else
        PLEs = ar.ple;
        PLEs.donePars = {};
        PLEs = loadSingle( PLEs, directory, filterProfile, appendProfiles);
    end
    
    PLEs.p_labels = ar.pLabel;
    PLEs.fighandel = zeros(size(ar.pLabel));
    
    if nargout == 0
        ar.ple = PLEs;
    else
        vargout = PLEs;
    end
end

function success = tryArLoad( PLEs, directory )
    global ar;
    
    success = 0; strict = 0;
    wsLoc = [ directory '/../../output.mat' ];
    if ( isfield( ar, 'pLabel' ) )
        if strict
                strictList = setdiff(ar.pLabel, PLEs.p_labels);
                nonStrict = setdiff(PLEs.p_labels, ar.pLabel);
        else
            strictList = [];
            nonStrict = setdiff(PLEs.p_labels, ar.pLabel);
        end
    else
        strictList = 1;
        nonStrict = 1;
    end
    if ( ( ( numel( strictList ) + numel( nonStrict ) > 0 ) )|| ( isempty( ar ) ) )
        fprintf( 'Incorrect or no model, loading struct ...' );
        if ( exist( wsLoc, 'file' ) )
            tmp = load( wsLoc, 'ar' );
            ar = tmp.ar;
            
            if ( ( ( numel( setdiff(ar.pLabel, PLEs.p_labels) ) + numel( setdiff(PLEs.p_labels, ar.pLabel) ) ) > 0 ) || ( isempty( ar ) ) )
                success = 1;
                disp( '[ OK ]\n' );
            end
        else
            error( 'Could not find model file' );
        end
    else
        success = 1;
    end
end

% Check if something is a ple struct
function out = isPLEStruct(ple)
    if ( isfield( ple, 'chi2s' ) && isfield( ple, 'ps' ) && isfield( ple, 'psinitstep' ) && isfield( ple, 'dchi2_point' ) )
        out = 1;
        return;
    else
        out = 0;
        return;
    end
end

% Directory contains either a directory to be loaded, or an ar struct to be
% loaded from.
function PLEs = loadSingle( PLEs, directory, doFilterProfile, doAppendProfiles)   
    global ar;
    strict = 0;
    
    if ischar(directory)
        pleLoc = [ directory '/PLE/results.mat' ];

        % Check whether a PLE is stored for this profile
        if ( ~exist( pleLoc, 'file' ) )
            disp( 'Empty dir' );
            return
        end
        tmp = load( pleLoc, 'ple' );
        curPLE = tmp.ple;

        % Check whether it is the correct model, otherwise reload the thing
        if ( tryArLoad( curPLE, directory ) == 0 )
            PLEs = ar.ple;
        end
    else
        if ( isfield( directory, 'ple' ) )
            directory = directory.ple;
        end
        if ~isPLEStruct( directory )
            help arMergeProfiles;
            error( 'Directory argument must contain either a directory that has PLEs, a PLE struct, an ar struct containing a PLE struct.' );
        end
        
        curPLE = directory;
    end
       
    % Start copying the profile
    filled = find( ~cellfun(@isempty, curPLE.chi2s) );
    notfound = [];
    
    if ~strict
        for b = 1 : numel(filled)
            idx = find( strcmp( ar.pLabel, curPLE.p_labels{filled(b)} ) );
            if ~isempty(idx)
                loc(b) = idx;
            else
                notfound = union( notfound, filled(b) );
            end
            
            if ~isempty( notfound )
                nf_str = sprintf('%s ', curPLE.p_labels{notfound});
                warning( 'Profiles present in profile struct that are not in the currently loaded model:\n%s\n', nf_str );
            
                % Remove profiles that have not been found in the target model
                % from the inclusion list
                filled(filled==notfound) = [];
            end
        end
    else
        loc = filled;
    end
    
    % Remap parameters
    [available, idxInCurPLE] = ismember( ar.pLabel, curPLE.p_labels );
    
    for b = 1 : length( filled )
        chi2s = curPLE.chi2s{filled(b)} + 0;
        if ( doFilterProfile )
            N = filterProfile(chi2s);
        else
            N = 1 : numel( chi2s );
        end
        fprintf( 'Found profile for %s\n', curPLE.p_labels{filled(b)} );
        
        appending = false;
        if sum(strcmp(PLEs.donePars, curPLE.p_labels{filled(b)})) ~= 0 && doAppendProfiles
            appending = true;
            fprintf('Combining profiles...\n')
        end
        PLEs.donePars{end + 1} = curPLE.p_labels{filled(b)};
        
        if appending
            PLEs = appendVectors(  PLEs, curPLE, {'chi2s', 'chi2sinit', 'chi2spriors', 'chi2spriorsAll', 'chi2sviolations'}, filled(b), loc(b), N );
            PLEs = appendMatrices( PLEs, curPLE, {'psinit', 'ps', 'gradient', 'psinitstep'}, filled(b), loc(b), N, available, idxInCurPLE );
            PLEs = copyScalars(  PLEs, curPLE, {'estimatetime', 'fittime', 'timing', 'conf_lb', 'conf_ub', 'conf_lb_point', 'conf_ub_point', 'IDstatus', 'p'}, filled(b), loc(b) );
            PLEs = copySingle(   PLEs, curPLE, {'breakon_point', 'dchi2', 'chi2_strID_ratio', 'initstep_fkt', 'minstepsize', 'breakonlb', 'breakonub', 'maxstepsize', ...
                'plot_point', 'plot_simu', 'dist_thres', 'grad_thres', 'dchi2_point', 'merit', 'alpha_level', 'ylabel', ...
                'integrate_fkt', 'fit_fkt', 'setoptim_fkt', 'merit_fkt', 'optimset_tol', 'allowbetteroptimum', 'savePath', 'relchi2stepincrease'} );
            PLEs.plot_point = curPLE.plot_point;
        else
            PLEs = copyVectors(  PLEs, curPLE, {'chi2s', 'chi2sinit', 'chi2spriors', 'chi2spriorsAll', 'chi2sviolations'}, filled(b), loc(b), N );
            PLEs = copyMatrices( PLEs, curPLE, {'psinit', 'ps', 'gradient', 'psinitstep'}, filled(b), loc(b), N, available, idxInCurPLE );
            PLEs = copyScalars(  PLEs, curPLE, {'estimatetime', 'fittime', 'timing', 'conf_lb', 'conf_ub', 'conf_lb_point', 'conf_ub_point', 'IDstatus',  'p'}, filled(b), loc(b) );
            PLEs = copySingle(   PLEs, curPLE, {'breakon_point', 'dchi2', 'chi2_strID_ratio', 'initstep_fkt', 'minstepsize', 'breakonlb', 'breakonub', 'maxstepsize', ...
                'plot_point', 'plot_simu', 'dist_thres', 'grad_thres', 'dchi2_point', 'merit', 'alpha_level', 'ylabel', ...
                'integrate_fkt', 'fit_fkt', 'setoptim_fkt', 'merit_fkt', 'optimset_tol', 'allowbetteroptimum', 'savePath', 'relchi2stepincrease'} );
            PLEs.plot_point = curPLE.plot_point;
        end
    end
end

function PLEs = copyVectors( PLEs, curPLE, names, ID, ID_target, N )
    for j = 1 : length( names )
        if ~isfield( PLEs, names{j} )
            PLEs.(names{j}) = {};
        end
        PLEs.(names{j}){ID_target} = curPLE.(names{j}){ID}(N) + 0;
    end
end

function PLEs = copyMatrices( PLEs, curPLE, names, ID, ID_target, N, available, idxInCurPLE )
    for j = 1 : length( names )
        if ~isfield( PLEs, names{j} )
            PLEs.(names{j}){ID_target} = [];
        end
        PLEs.(names{j}){ID_target}(:, available) = curPLE.(names{j}){ID}(N,idxInCurPLE(available)) + 0;
    end
end

function PLEs = copyScalars( PLEs, curPLE, names, ID, ID_target )
    for j = 1 : length( names )
        PLEs.(names{j})(ID_target) = curPLE.(names{j})(ID) + 0;
    end
end

function PLEs = copySingle( PLEs, curPLE, names )
    for j = 1 : length( names )
        PLEs.(names{j}) = curPLE.(names{j});
    end
end

function PLEs = appendVectors( PLEs, curPLE, names, ID, ID_target, N )
    for j = 1 : length( names )
        if ~isfield( PLEs, names{j} )
            PLEs.(names{j}){ID_target} = zeros(1, numel(available));
        end
        PLEs.(names{j}){ID_target}(end+1:end+length(N)) = curPLE.(names{j}){ID}(N) + 0;
    end
end

function PLEs = appendMatrices( PLEs, curPLE, names, ID, ID_target, N, available, idxInCurPLE )
    for j = 1 : length( names )
        if ~isfield( PLEs, names{j} )
            PLEs.(names{j}){ID_target} = [];
        end
        PLEs.(names{j}){ID_target}(end+1:end+length(N),available) = curPLE.(names{j}){ID}(N,idxInCurPLE(available)) + 0;
    end
end

function N = filterProfile(chi2s)
    done = 0;
    fprintf( '%d =>', numel(chi2s) );
    remaining = 1:numel(chi2s);
    remaining( isnan(chi2s) ) = [];
    while( ~done )
        [N, cull] = findPeaks(chi2s(remaining));
        if ( ~isempty(cull) )
            remaining(cull) = [];
        else
            done = 1;
        end
    end
    
    N = remaining;
    fprintf( ' %d\n', numel(N) );
end

function [N, cull] = findPeaks( signal )
    dc = diff( signal );
    sn = sign(dc(1:end-1))>sign(dc(2:end));
    [N,I] = sort( dc, 'descend' );
    oo5 = N(floor(length(dc)*0.2));
    
    % Remove the ones with biggest derivative
    med = median(dc);
    md  = mad(dc, 1);
    cull = intersect(find(sn)+1, find(dc>oo5)+1);

    mask = ones(numel(signal), 1); cull = [];
    mask( cull ) = 0;   
    N = find(mask>0);
end

