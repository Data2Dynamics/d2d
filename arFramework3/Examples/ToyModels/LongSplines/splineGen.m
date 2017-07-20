% This function generates an input for a monotonic spline of a desired
% length. Note that it is not optimal in the sense that it pads the spline
% up to a certain number of points such that it fits into monospline10s.
% One could edit the final result by hand to get a better match with
% whatever one is fitting.

function spline = splineGen(t_in)    
    tlen = numel(t_in);
    tlen = tlen - 8;                    % First spline can support 8 points
    tlen = ceil( tlen / 7 ) * 7 + 8;    % Subsequent splines 7
    t = nan(tlen, 1);
    t(1:numel(t_in)) = t_in;

    knot_pattern = {};
    t_pattern = [];
    epsy = 0.1;                         % Distance between same knots
    ct = 1;                             % Keeps track of current time point
    knotID = 1;                         % Keeps track of position in current spline
    currentKnot = 1;                    % Current knot parameter
    q = 1;                              % Keeps track of pattern index
    while (ct <= numel(t))
        cKnotString = sprintf('knot%d', currentKnot);

        % Start of a spline
        if ( knotID == 1 )
            if ( q > 1 )
                cKnotString = sprintf( '%s - %s', cKnotString, knot_pattern{q-1} );
            end
            [knot_pattern{q}, t_pattern(q)] = addKnot( q, knot_pattern, t_pattern, cKnotString, t(ct) - epsy );
            q = q + 1;
        end

        % Knots
        [knot_pattern{q}, t_pattern(q)] = addKnot( q, knot_pattern, t_pattern, cKnotString, t(ct) );
        knotID = knotID + 1;
        currentKnot = currentKnot + 1;
        q = q + 1;

        % End of a spline
        if ( knotID > 8 )
            knotID = 1;
            [knot_pattern{q}, t_pattern(q)] = addKnot( q, knot_pattern, t_pattern, cKnotString, t(ct) + epsy );
            q = q + 1;
        else
            ct = ct + 1;
        end
    end

    t_pattern = t_pattern(1:end-2);
    knot_pattern = knot_pattern(1:end-2);
    
    % Build the string
    c = 1;
    for q = 1 : 10 : numel( t_pattern )
        splineFunc{c} = 'monospline10(t';
        for r = 1 : 10
            splineFunc{c} = sprintf('%s, %d, %s', splineFunc{c}, t_pattern(q+r-1), knot_pattern{q+r-1});
        end
        splineFunc{c} = sprintf('%s)', splineFunc{c});
        c = c + 1;
    end

    % Add an offset, such that the spline coefficients don't have to go
    % negative. This is desirable since fitting in log space is better!
    spline = sprintf('%s + ', splineFunc{:});
    spline = sprintf('%s - spline_offset', spline(1:end-2));
end

function [knot, t] = addKnot( q, knot_pattern, t_pattern, knot, t )
    if ( isnan(t) )
        t = t_pattern(q-1) + 0.01;
        knot = 'spline_dummy';
    end
end