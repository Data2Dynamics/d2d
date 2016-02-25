% arInitialObservationEstimate( frac, datasets, additional arguments )
%
% Function to prevent solutions "ignoring" specific datasets when fitting 
% with error model. Computes and sets an upper bound for the estimated
% error. Optionally roughly estimates observation parameters such as scale 
% and offset for models with errors of the form scale * X + offset.
%
% Parameters:
%   frac                - Percentage of unexplained variance which is O.K.
%                         Default: 0.25
%   datasets            - List of data sets to estimate error (and optionally observational parameters for)
%   
% Flags (supply as string)
%   estimateObsPars     - Estimate observational parameters
%   symbolic            - Use this if your observational parameter names do not start with scale or offset. Much slower however.
%   verbose             - Show the determined values

function arInitialObservationEstimate( varargin )

    global ar;
    
    try
        arSimu(false,false,true);
    catch
    end

    estimateVariance = 1;   % Do we need to estimate variances?
    frac = 0.25;            % Percentage of unexplained variance which is considered O.K.
    
    if (~isempty(varargin)) && ( isnumeric( varargin{1} ) )
        frac        = varargin{1};
        varargin    = varargin(2:end);
        arFprintf(2, 'Fraction set to %.3f\n', frac );
    end
    
    if (~isempty(varargin)) && ( isnumeric( varargin{1} ) )
        datasets        = varargin{1};
        varargin        = varargin(2:end);
        datasetsString  = sprintf( '%d ', datasets );
        arFprintf(2, 'Datasets set to %s\n', datasetsString );
    end
    
    opts = argSwitch( {'estimateobspars', 'symbolic', 'verbose'}, varargin{1:end} );
    
    m = 1;   
    if ~exist( 'datasets', 'var' )
        datasets = 1 : length( ar.model(m).data );
    end

    if ( opts.symbolic )
        vX = cellsym(ar.model(m).x);
        vZ = cellsym(ar.model(m).z);
        vA = cellsym(union(ar.model(m).x, ar.model(m).z));
    end
    
    stdLists        = struct;
    scaleLists      = struct;
    scaleSimLists   = struct;
    offsetLists     = struct;
    offsetSimLists  = struct;
    parNames        = {};
    scaleNames      = {};
    offsetNames     = {};
    for jdn = 1 : length( datasets )
        jd = datasets(jdn);

        cond  = ar.model(m).data(jd).cLink;
        tPts  = ar.model(m).data(jd).tLinkExp;
        if tPts < 2
            tSelect = 1 : length( ar.model(m).condition(cond).tExp );
        else
            tSelect = tPts;
        end

        zVars = ar.model(m).condition( cond ).zExpSimu(tSelect, :);
        xVars = ar.model(m).condition( cond ).xExpSimu(tSelect, :);

        stdev = cellsym( ar.model(m).data(jd).fystd );
        if ( opts.symbolic )
            stdev = cellsubs(stdev, vX, zeros(size(vX)));
            stdev = cellsubs(stdev, vZ, zeros(size(vZ)));    
        end
        
        for jy = 1 : length( ar.model(m).data(jd).fy )
            if ( opts.estimateobspars )
                if ( opts.symbolic )
                    offset = subs(ar.model(m).data(jd).fy{jy}, vA, zeros(size(vA)));
                    scale  = subs(ar.model(m).data(jd).fy{jy}, vA, ones(size(vA)))-offset;
    
                    if ( scale == 0 )
                        scale = [];
                    end
                    if ( offset == 0 )
                        offset = [];
                    end
                else
                    [s,t] = regexp(ar.model(m).data(jd).fy{jy}, 'offset_(\w*)');
                    offset = ar.model(m).data(jd).fy{jy}(s:t);
                    [s,t] = regexp(ar.model(m).data(jd).fy{jy}, 'scale_(\w*)');
                    scale = ar.model(m).data(jd).fy{jy}(s:t);
                end
                                
                % We have a scale to estimate
                if ~isempty( scale )
                    if ( ar.model(m).data(jd).logfitting(jy) )
                        scaleLists      = addToList( scaleLists, scale, 10.^ar.model(m).data(jd).yExp(:,jy) );
                        scaleSimLists   = addToList( scaleSimLists, scale, 10.^ar.model(m).data(jd).yExpSimu(:,jy) );
                    else
                        scaleLists      = addToList( scaleLists, scale, ar.model(m).data(jd).yExp(:,jy) );
                        scaleSimLists   = addToList( scaleSimLists, scale, ar.model(m).data(jd).yExpSimu(:,jy) );
                    end
                    scaleNames      = union( scaleNames, char(scale) );
                    
                    % We have an offset to estimate
                    if ~isempty( offset )
                        if ( ar.model(m).data(jd).logfitting(jy) )
                            offsetLists = addToList( offsetLists, offset, 10.^ar.model(m).data(jd).yExp(:,jy) );
                            offsetSimLists = addToList( offsetSimLists, offset, 10.^ar.model(m).data(jd).yExpSimu(:,jy) );
                        else
                            offsetLists = addToList( offsetLists, offset, ar.model(m).data(jd).yExp(:,jy) );
                            offsetSimLists = addToList( offsetSimLists, offset, ar.model(m).data(jd).yExpSimu(:,jy) );
                        end
                        
                        offsetNames = union( offsetNames, char(offset) );
                    end
                end
            end
            
            % If we can compute a meaningful variance from this, subtract
            % the mean and add it to the list corresponding to that
            % parameter
            if ( estimateVariance )
                cvar = nanvar(ar.model(m).data(jd).yExp(:,jy));
                if ~isnan(cvar)
                    parNames = union( parNames, char(stdev{jy}) );
                    
                    stdLists = addToList( stdLists, stdev{jy}, ar.model(m).data(jd).yExp(:,jy) - mean(ar.model(m).data(jd).yExp(:,jy)) );
                end
            end
        end
    end
    
    for a = 1 : length( parNames )
        sdMax = sqrt( frac * nanvar(stdLists.(hashstring(parNames{a}))) );
        pID = find( ismember(ar.pLabel, parNames{a} ) );
        qLog10 = ar.qLog10(pID);
        if ( qLog10 )
            sdMax = log10(sdMax);
        end
        ar.ub(pID) = sdMax;
        ar.p(pID) = sdMax;
        
        if ( opts.verbose )
            arFprintf( 2, '%s: %d\n', parNames{a}, ar.p(pID) );
        end        
    end
    
    if ( opts.estimateobspars )
        % Estimate the scales
        for a = 1 : length( scaleNames )
            chash = hashstring(scaleNames{a});
            maxValExp   = nanmax(scaleLists.(chash));
            minValExp   = nanmin(scaleLists.(chash));
            
            maxValSim   = nanmax(scaleSimLists.(chash));
            minValSim   = nanmin(scaleSimLists.(chash));
            ratio       = (maxValExp - minValExp) / (maxValSim - 0 );
                       
            % If it's a valid correction, do it
            if (~isnan(ratio))
                pID = find( ismember(ar.pLabel, scaleNames{a} ) );
                if ( ar.qLog10(pID) )
                    ar.p(pID) = ar.p(pID) + log10(ratio);
                else
                    ar.p(pID) = ar.p(pID)*ratio;
                end
                
                if ( opts.verbose )
                    arFprintf( 2, '%s: %d\n', scaleNames{a}, ar.p(pID) );
                end 
            end
        end
                
        % Estimate the offsets
        for a = 1 : length( offsetNames )
            chash = hashstring(offsetNames{a});
            if numel( offsetLists.(chash) ) > 1
                minOffsetExp = nanmin(offsetLists.(chash));
                minOffsetSim = nanmin(offsetSimLists.(chash));
            else
                minOffsetExp = 10^-5;
                minOffsetSim = 10^-5;
            end
            
            pID = find( ismember(ar.pLabel, offsetNames{a} ) );
            
            if ( ar.qLog10(pID) )
                ar.p(pID) = ar.p(pID) - log10(minOffsetSim) + log10(minOffsetExp);
            else
                ar.p(pID) = minOffsetExp;
            end
            
            if ( opts.verbose )
               arFprintf( 2, '%s: %d\n', offsetNames{a}, ar.p(pID) );
            end 
        end     
    end
    
    for pID = 1 : length( ar.p )
        ar.p(pID) = clamp(ar.p(pID), ar.lb(pID), ar.ub(pID));
    end
    
end

function f = clamp(f, mi, ma)
    f = max( [ min( [f, ma] ), mi ] );
end

function strct = addToList( strct, listfield, values )
    curField = hashstring(listfield);
    if ( ~isfield( strct, curField ) )
        strct.(curField) = values;
    else
        strct.(curField) = [ strct.(curField) ; values ];
    end
end

function cellarray = cellsym( list )
    cellarray = cell(size(list));
    for a = 1 : length( list )
        cellarray{a} = sym(list{a});
    end
end

function list = cellsubs( list, old, new )
    for a = 1 : length( list )
        list{a} = subs(list{a}, old, new);
    end
end

function hash = hashstring(str)
    algs = {'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
    checksum = java.security.MessageDigest.getInstance(algs{2});
    str = char(str);
    checksum.update(uint8(str(:)));
    h = typecast(checksum.digest,'uint8');
    hash = dec2hex(h)';
    hash = ['F', hash(:)'];
    clear checksum;
end

function [opts] = argSwitch( switches, varargin )

    for a = 1 : length(switches)
        opts.(switches{a}) = 0;
    end

    for a = 1 : length( varargin )
        if ( max( strcmpi( varargin{a}, switches ) ) == 0 )
            error( 'Invalid switch argument was provided %s', varargin{a} );
        end
        opts.(lower(varargin{a})) = 1;
    end    
end