% [used, m] = arFindParameterUse(idx, [verbose] )
%
% Find data sets and conditions in which a certain parameter is used
%
%    idx            index of parameter
%    verbose  [0]   boolean, specifies if output is displayed
%
%    used           boolean, if parameter ar.p(idx) is used at all
%    m              struct, with fields
%                           ci: list with all conditions containing the parameter ar.p(idx)
%                           di contains list with all data containing the parameter ar.p(idx)

function [used, m] = arFindParameterUse( idx, verbose )

    if ( nargin < 2 )
        verbose = 0;
    end

    global ar;
    if ( nargout > 0 )
        used = 0;
    end
    for jm = 1 : numel( ar.model )
        if ( nargout > 1 )
            m(jm) = struct;
            m(jm).ci = [];
            m(jm).di = [];
        end
        for jc = 1 : numel( ar.model(jm).condition )
            if ( ~isempty( ar.model(jm).condition(jc).pLink ) && ar.model(jm).condition(jc).pLink(idx) )
                if ( verbose )
                    fprintf( 'Used in condition %d\n', jc );
                end
                if ( nargout > 1 )
                    m(jm).ci = union( m(jm).ci, jc );
                end
                if ( nargout > 0 )
                    used = 1;
                end
            end
        end
        for jd = 1 : numel( ar.model(jm).data )
            if ( ~isempty( ar.model(jm).data(jd).pLink ) && ar.model(jm).data(jd).pLink(idx) )
                if ( verbose )
                    fprintf( 'Used in data %d: %s\n', jd, ar.model(jm).data(jd).name );
                end
                if ( nargout > 1 )
                    m(jm).di = union( m(jm).di, jd );
                end
                if ( nargout > 0 )
                    used = 1;
                end
            end
        end
    end