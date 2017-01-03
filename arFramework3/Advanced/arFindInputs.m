% arFindInputs
%
% Usage: 
%    arFindInputs;
%
% Finds keywords step1 and step2 in the input functions and adds the time
% points they correspond to as events in the ar structure.
%
% Currently, only events with either numeric time points, parameters or 
% condition values are supported. Events that do not correspond to 
% this format raise a warning but are ignored.
%
% Event locations are stored in ar.model(#).condition(#).tEvents
% Please note that a relink of the model is required after running this 
% command to be able to make use of the event system (arLink).
% An additional recompile is *not* required.
%
% Events can be enabled or disabled globally by setting ar.config.useEvents
% to 0 or 1; or locally by setting ar.model(#).condition(#).qEvents to 0 
% after linking. Note that if an event is set to zero manually, it will not
% automatically be re-enabled upon the next linkage. To re-enable the
% events for a specific condition, manually set ar.model(#).condition(#).qEvents
% to 1.

function arFindInputs( verbose )

    global ar;
    global arOutputLevel;

    if ( nargin < 1 )
        verbose = 0;
    end

    logCall( 'arFindInputs', verbose );

    % Sort parameter labels by size to avoid replacing ones which are a shorter
    % substring of longer ones first.
    [~, I] = sort( cellfun( @length, ar.pLabel ), 'descend' );

    % Find out how long to make the waitbar
    for m = 1 : length( ar.model )
        allData = length( ar.model(m).condition );
    end

    if ( arOutputLevel > 1 )
        h = waitbar(0);
    end
    lastParse = ''; lastA = 0; totalEvents = 0;
    for m = 1 : length( ar.model )
        for a = 1 : length( ar.model(m).condition )
            if ( arOutputLevel > 1 )
                waitbar(a/allData, h, sprintf( 'Processing step inputs [%d/%d] %s', a, allData, lastParse ) );
            end

            events = [];
            stepLocations = {};
            % Find the event time points
            for b = 1 : length( ar.model(m).condition(a).fu );
                step1 = findstr(ar.model(m).condition(a).fu{b}, 'step1'); %#ok
                step2 = findstr(ar.model(m).condition(a).fu{b}, 'step2'); %#ok

                for c = 1 : length( step1 )
                    ar.model(m).condition(a).fu{b}(step1(c):end);
                    chk = strsplit(ar.model(m).condition(a).fu{b}(step1(c):end),',');
                    stepLocations{end+1} = chk{3}; %#ok
                end
                for c = 1 : length( step2 )
                    chk = strsplit(ar.model(m).condition(a).fu{b}(step2(c):end),',');
                    stepLocations{end+1} = chk{3}; %#ok
                    stepLocations{end+1} = chk{5}; %#ok
                end        
            end

            % Transform the parameters that are defined in log space
            parVals = ar.p;
            parVals( ar.qLog10 > 0 ) = 10.^parVals( ar.qLog10 > 0 );

            for b = 1 : length( stepLocations )
                l = str2num(stepLocations{b}); %#ok
                
                % Not a numeric value, check whether parameter or condition 
                % substitution leads to a valid expression
                if isempty( l )
                    l       = stepLocations{b};
                    lold    = l;

                    % Is it in the condition variable list? ==> Replace it
                    [~,I2] = sort(cellfun(@length,ar.model(m).condition(a).pold), 'descend');
                    pold = ar.model(m).condition(a).pold(I2);
                    pnew = ar.model(m).condition(a).fp(I2).';
                    for c = 1 : length( pold )
                        l = strrep( l, pold{c}, pnew{c} );
                    end
                    
                    % Parameter values
                    for c = 1 : length( ar.pLabel )
                        l = strrep( l, ar.pLabel{I(c)}, num2str(parVals(I(c))) );
                    end
                    
                    l = str2num(l); %#ok
                    lastParse = sprintf( '%s => %s in condition %d', strrep( lold, '_', '\_' ), strrep( num2str(l), '_', '\_' ), a );
                end

                try
                    events(end+1) = l; %#ok
                    if ( verbose )
                        arFprintf(1, l);
                    end
                catch
                    arFprintf( 1, 'Failure to parse input %s for model %d and data %d', stepLocations{b}, m, a );
                end
            end

            totalEvents = totalEvents + length( events );

            if isfield( ar.model(m).condition(a), 'tEvents' )
                ar.model(m).condition(a).tEvents = union( ar.model(m).condition(a).tEvents, events );
            else
                ar.model(m).condition(a).tEvents = union( [], events );
            end
        end
        lastA = lastA + length( ar.model(m).condition );
    end

    arFprintf( 1, '%d input events assigned!\n', totalEvents );
    if ( arOutputLevel > 1 )
        close(h);
    end

    ar.config.useEvents = 1;
    arLink(true);
    
    % Invalidate cache so simulations do not get skipped
    arCheckCache(1);
end


function logCall( fn, varargin )
    global ar;
    
    if ~isfield(ar, 'eventLog')
        ar.eventLog = {};
    end
    call = [fn '('];
    if length( varargin ) > 0 %#ok
        call = sprintf('%s%s', call, mat2str(varargin{1}) );
    end
    for a = 2 : length( varargin )
        call = sprintf('%s, %s', call, mat2str(varargin{a}) );
    end
    call = [call ')'];
    ar.eventLog{length(ar.eventLog)+1} = call;
end