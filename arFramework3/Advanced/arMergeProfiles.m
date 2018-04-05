%
% arMergeProfiles( directory, forceLoad )
%
%   Merge profiles from various runs
%
%   Arguments:
%     - Directory
%         Should point to results directory. Typically, the profiles will be
%         subdirectories in this directory; each with a subdirectory /PLE/
%     - forceLoad
%         Force reload of the model
%     - filterProfile
%         Smoothen profile? (default: 0 (off))
%
%  If no model is loaded, an incompatible model is loaded or forceLoad is
%  set to one, arMergeProfiles will attempt to load a model. For this, it
%  will look one directory up from Directory for a mat file named
%  arStruct.mat
% 
function vargout = arMergeProfiles( directory, forceLoad, filterProfile )

    global ar;
    s = dir(directory);

    if ( ~exist( 'forceLoad', 'var' ) )
        forceLoad = 0;
    end
    if ( ~exist( 'filterProfile', 'var' ) )
        filterProfile = 0;
    end
    
    if isempty( ar ) || forceLoad
        load( [ directory '/../arStruct.mat' ] );
    end
    
    % Assemble ple struct
    PLEs = ar.ple;
    for a = 1 : length( s )
        if ( ~strcmp( s(a).name, '.' ) && ~strcmp( s(a).name, '..' ) )
            path = [ directory '/' s(a).name ];
            if ( isdir( path ) )
                PLEs = loadSingle( PLEs, path, filterProfile );
            end
        end
    end
    
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

function PLEs = loadSingle( PLEs, directory, doFilterProfile )   
    global ar;
    strict = 0;
    
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
       
    % Start copying the profile
    filled = find( ~cellfun(@isempty, curPLE.chi2s) );
    
    if ~strict
        for b = 1 : numel(filled)
            loc(b) = find( strcmp( ar.pLabel, curPLE.p_labels{filled(b)} ) );
            mapping(b) = find( strcmp( ar.pLabel, curPLE.p_labels{filled(b)} ) );
        end
    else
        loc = filled;
    end
    
    for b = 1 : length( filled )
        chi2s = curPLE.chi2s{filled(b)} + 0;
        if ( doFilterProfile )
            N = filterProfile(chi2s);
        else
            N = 1 : numel( chi2s );
        end
        fprintf( 'Found profile for %s\n', curPLE.p_labels{filled(b)} );
        PLEs = copyVectors(  PLEs, curPLE, {'chi2s', 'chi2sinit', 'chi2spriors', 'chi2spriorsAll', 'chi2sviolations'}, filled(b), loc(b), N );
        PLEs = copyMatrices( PLEs, curPLE, {'psinit', 'ps', 'gradient', 'psinitstep'}, filled(b), loc(b), N );
        PLEs = copyScalars(  PLEs, curPLE, {'estimatetime', 'fittime', 'timing', 'conf_lb', 'conf_ub', 'conf_lb_point', 'conf_ub_point'}, filled(b), loc(b) );
        PLEs = copySingle(   PLEs, curPLE, {'plot_point', 'plot_simu', 'dist_thres', 'grad_thres', 'dchi2_point', 'merit', 'alpha_level', 'p', 'p_labels', 'ylabel'} );
        PLEs.plot_point = curPLE.plot_point;
    end
end

function PLEs = copyVectors( PLEs, curPLE, names, ID, ID_target, N )
    for j = 1 : length( names )
        try
            PLEs.(names{j}){ID_target} = curPLE.(names{j}){ID}(N) + 0;
        catch
        end
    end
end

function PLEs = copyMatrices( PLEs, curPLE, names, ID, ID_target, N )
    for j = 1 : length( names )
        try
            PLEs.(names{j}){ID_target} = curPLE.(names{j}){ID}(N,:) + 0;
        catch
        end
    end
end

function PLEs = copyScalars( PLEs, curPLE, names, ID, ID_target )
    for j = 1 : length( names )
        try
            PLEs.(names{j})(ID_target) = curPLE.(names{j})(ID) + 0;
        catch
        end
    end
end

function PLEs = copySingle( PLEs, curPLE, names )
    for j = 1 : length( names )
        PLEs.(names{j}) = curPLE.(names{j});
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

