% Work in progress ...
% Full description will follow after further testing
% 
% function [out, outL, L2, t] = arParameterScan(parameter, inRange, predictions, varargin)

function [out, outL, L2, t] = arParameterScan(parameter, inRange, predictions, varargin)

    global arScanParameters;
    global ar;
    
    switches = { 'log', 'model', 'condition', 't', 'steps' };
    extraArgs = [ 1, 1, 1, 1, 1 ];
    description = { ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} ...
    {'', ''} };
    [opts, varargin] = argSwitch( switches, extraArgs, description, 0, 'softmatching', varargin );    
    
    if ( opts.log )
        scanInLog = opts.log_args;
    else
        scanInLog = 1;
    end
    if ( opts.model )
        m = opts.model_args;
    else
        m = 1;
    end    
    if ( opts.condition )
        c = opts.condition_args;
    else
        c = 1;
    end
    if ( opts.t )
        if ( opts.t_args )
            [tFine, iFine] = intersect( ar.model(m).condition(c).tFine, opts.t_args );
            [tExp, iExp]  = intersect( ar.model(m).condition(c).tExp,  opts.t_args );
            t = union( tFine, tExp );
            [~, iA] = intersect( t, tFine );
            [~, iB] = intersect( t, tExp );
            
            if ( isempty( t ) )
                error( 'Specified time point(s) do not exist in tFine' );
            end
        end
    else
        tExp  = ar.model(m).condition(c).tExp(end);
        tFine = ar.model(m).condition(c).tFine(end);
        if ( tExp > tFine )
            tFine = [];
            iFine = [];
            iA = [];
            iB = 1;
            iExp = numel(ar.model(m).condition(c).tExp);
        else
            tExp = [];
            iExp = [];
            iB = [];
            iA = 1;
            iFine = numel(ar.model(m).condition(c).tFine);
        end
    end
    
    t = union( tFine, tExp );
    if ( iscell( parameter ) )
        IDs = arFindPar( parameter, 'exact' );
    end
    if ( ischar( parameter ) )
        IDs = arFindPar( parameter, 'exact' );
    end
    if ( numel(IDs) > 2 )
        error( 'Function only supports 2 dimensional scans' );
    end
    if ( numel(IDs) == 2 )
        % Two variables
        if ( opts.steps )
            nsteps = opts.steps_args;
        else
            nsteps = 5;
        end        
        if ~iscell( inRange )
            range{1} = inRange;
            range{2} = inRange;
        end
    elseif ~iscell( inRange )
        % Single variable
        if ( opts.steps )
            nsteps = opts.steps_args;
        else
            nsteps = 20;
        end
        range{1} = inRange;
        range{2} = [1, 1];
    end
	if ( numel( IDs ) == 1 )
        IDs(2) = IDs(1);
    end
    
    lin = @(range)range(1) : max([1e-6, (range(2)-range(1))/(nsteps-1)]) : range(2);
    pLabel1 = ar.pLabel{IDs(1)};
    pLabel2 = ar.pLabel{IDs(2)};
    if ( scanInLog )
        % Log parameter ranges
        P1 = [ arGetPars( pLabel1, 1 ) + log10(range{1}(1)), arGetPars( pLabel1, 1 ) + log10(range{1}(2)) ];
        P2 = [ arGetPars( pLabel2, 1 ) + log10(range{2}(1)), arGetPars( pLabel2, 1 ) + log10(range{2}(2)) ];
        L1 = lin(P1);
        L2 = lin(P2);
    else
        % Linear parameter ranges
        P1 = [ arGetPars( pLabel1, 0 ) * range{1}(1), arGetPars( pLabel1, 0 ) * range{1}(2) ];
        P2 = [ arGetPars( pLabel2, 0 ) * range{2}(1), arGetPars( pLabel2, 0 ) * range{2}(2) ];
        L1 = lin(P1);
        L2 = lin(P2);
    end
    
    % Transform the range to whatever is used in the ar struct
    pv1 = trafo( ar.qLog10(IDs(1)), scanInLog, L1 );
    pv2 = trafo( ar.qLog10(IDs(1)), scanInLog, L2 );

    try
        changed = sum( abs( arScanParameters.pv1 - pv1 ) + abs( arScanParameters.pv2 - pv2 ) + abs( arScanParameters.m - m ) + abs( arScanParameters.c - c ) ) > 0;
    catch
        changed = 1;
    end
    
    % Calculations
    if ( changed )
        predictionsFine = cell( numel( L1 ), numel( L2 ) );
        predictionsExp  = cell( numel( L1 ), numel( L2 ) );
        for jx = 1 : numel( L1 )
            fprintf( 'Simulating %d/%d ...\n', jx, numel(L1) );
            for jy = 1 : numel( L2 )
                arPush('arParameterScan');
                try
                    % Do not switch these!
                    ar.p(IDs(2)) = pv2(jy);
                    ar.p(IDs(1)) = pv1(jx);
                    arSimu(false, true, true);
                    arSimu(false, false, true);
                    predictionsFine{jx,jy}  = [ar.model(m).condition(c).xFineSimu, ar.model(m).condition(c).zFineSimu, ar.model(m).condition(c).vFineSimu] + 0;
                    predictionsExp{jx,jy}   = [ar.model(m).condition(c).xExpSimu,  ar.model(m).condition(c).zExpSimu,  ar.model(m).condition(c).vExpSimu]  + 0;
                catch
                    predictionsFine{jx,jy}  = NaN * [ar.model(m).condition(c).xFineSimu, ar.model(m).condition(c).zFineSimu,ar.model(m).condition(c).vFineSimu] + 0;
                    predictionsExp{jx,jy}   = NaN * [ar.model(m).condition(c).xExpSimu,  ar.model(m).condition(c).zExpSimu, ar.model(m).condition(c).vExpSimu]  + 0;
                end
                arPop('silent');
            end
        end
        arScanParameters.pv1                = pv1;
        arScanParameters.pv2                = pv2;
        arScanParameters.m                  = m;
        arScanParameters.c                  = c;
        arScanParameters.predictionsFine    = predictionsFine;
        arScanParameters.predictionsExp     = predictionsExp;
    else
        disp( 'Using cached values' );
    end
    
    % Plotting
    modelVariables = [ar.model(m).x, ar.model(m).z, ar.model(m).v]; %#ok
    [~,I] = intersect( modelVariables, predictions );
    
    % Timepoints  
    % iFine => Indices of tFine to grab
    % iExp => Indices of iExp to grab
    % iA => location tFine result should be stored at
    % iB => location tExp result should be stored at
    out = NaN( numel(L1), numel(L2), numel(I), numel(t) );
    for jx = 1 : numel( L1 )
        for jy = 1 : numel( L2 )
            fine    = arScanParameters.predictionsFine{jx, jy};
            exp     = arScanParameters.predictionsExp{jx, jy};
            out(jx, jy, 1:numel(I), iA) = fine(iFine, I);
            out(jx, jy, 1:numel(I), iB) = exp(iExp, I);
        end
    end
    
    if ( ndims(IDs) > 2 )
        outL = {L1, L2};
    else
        outL = L1;
    end
end

function values = trafo( ar_log10, scanInLog, in )
    if ( ar_log10 )
        if ( scanInLog )
            values = in;
        else
            values = log10( in );
        end
    else
        if ( scanInLog )
            values = 10.^in;
        else
            values = in;
        end
    end
end
    