% arDetectionLimit( m, d, LoD, [LL] )
%
% Specify a detection limit for a specific model m and dataset d.
%
%   m       -   Model index
%   d       -   Data index (you can find this index with arFindData)
%   LoD     -   Level of detection (vector of length ny, where ny is the
%               number of observables in this dataset). Set to -inf 
%               for observables that are not hampered by a detection limit.
%   LL      -   Lower limit of data censoring [-inf])
%
% Detection limits are implemented by including a term corresponding to the
% likelihood integrated up to the detection limit for those values.
%
% Example:
%   arDetectionLimit( 1,1,[-inf,5,-inf])
%   if observable has a lower detection limit of 5
%

function arDetectionLimit( m, d, LoD, LL )
    global ar;
    
    if ( nargin < 4 )
        LL = -inf(size(LoD));
    end
    ar.model(m).data(d).resfunction.type                = 'DetectionLimit';
    ar.model(m).data(d).resfunction.active              = 1;
    ar.model(m).data(d).resfunction.LoD                 = LoD;
    ar.model(m).data(d).resfunction.fres_fun            = @(yexp, y, ystd, fiterrors_correction_factor)fres_LoD(yexp, y, ystd, fiterrors_correction_factor, LoD, LL);
    ar.model(m).data(d).resfunction.fres_error_fun      = @(yexp, ystd, add_c)fres_error_LoD(yexp, ystd, add_c, LoD, LL);
    ar.model(m).data(d).resfunction.fsres_fun           = @(yexp, y, sy, ystd, systd, res, reserr, add_c, fiterrors_correction_factor)fsres_LoD(yexp, y, sy, ystd, systd, res, reserr, add_c, fiterrors_correction_factor, LoD, LL);
    ar.model(m).data(d).resfunction.fsres_error_fun     = @(yexp, y, sy, ystd, systd, res, reserr, add_c, fiterrors_correction_factor, sres)fsres_error_LoD(yexp, y, sy, ystd, systd, res, reserr, add_c, fiterrors_correction_factor, sres, LoD, LL);
end

% Least squares with detection limit handling
%
% - log likelihood = N_aboveLoD log(sigma) + 0.5 * sum( ( (d-f)/sig )^2 ) - sum( log( Phi( (LoD-y)/sig ) - Phi( (LL-y)/sig ) )
% We work with -2 log likelihood, hence the factor 2.
%
function [res,chi2] = fres_LoD(yexp, y, ystd, fiterrors_correction_factor, LoD, LL)
    LoD = repmat( LoD, size(yexp,1), 1 );
    LL  = repmat( LL, size(yexp,1), 1 );
    
    idx = yexp <= LoD;
    yidx = y(idx);
    ystdidx = ystd(idx);
    LoDidx = LoD(idx);
    LLidx = LL(idx);
    
    res = NaN( size( y ) );
    
    % Call normal fres for the ones above the detection limit
    % The eps(1) is for numerical reasons. It introduces a negligible bias
    % which helps prevent problems when the simulations are very far from
    % the detection limit.
    res(~idx) = fres( yexp(~idx), y(~idx), ystd(~idx), fiterrors_correction_factor );
    res(idx)  = sqrt( - 2 * log( eps(1) + Phi( ( LoDidx - yidx ) ./ ystdidx ) - Phi( ( LLidx - yidx ) ./ ystdidx ) ) );
    
    chi2 = sum( res.^2, 1 );
end

% Probability density integral of a standard normal from -infty to x
function norm_int = Phi(x)
    norm_int = 0.5 * ( 1 + erf( x / sqrt(2) ) );
end

% Standard normal
function stdn = phi(x)
    stdn = ( 1 / (sqrt(2*pi)) ) * exp( - (x.^2)/2 );
end

% Derivative of UDL w.r.t. p
function dudldp = diffUDL_p( y, sy, systd, ystd, LoD )
    inverseSD = repmat( 1 ./ ystd, 1, 1, size(sy,3) );
    UDLterm = repmat( ( LoD - y ) ./ ystd, 1, 1, size(sy,3) );

    dudldp = - inverseSD .* ( sy + systd .* UDLterm );    
end

% Derivative of UDL w.r.t. p
function dulldp = diffULL_p( y, sy, systd, ystd, LL )
    inverseSD = repmat( 1 ./ ystd, 1, 1, size(sy,3) );
    ULLterm = repmat( ( LL - y ) ./ ystd, 1, 1, size(sy,3) );
    
    dulldp = - inverseSD .* ( sy + systd .* ULLterm );    
end

% Via the chain rule we can obtain:
%   - ( phi(UDL)*diff(UDL, p) - phi(ULL)*diff(ULL, p) ) / ( res * ( Phi(UDL) - Phi(ULL) ) )
% With res being defined as before (including the square root).
%
% The derivatives of UDL and ULL can potentially contain a derivative with respect to sigma.
function sres = fsres_LoD(yexp, y, sy, ystd, systd, res, reserr, add_c, fiterrors_correction_factor, LoD, LL)
    LoD     = repmat( LoD, size(yexp,1), 1 );
    LL      = repmat( LL, size(yexp,1), 1 );
    idx     = yexp <= LoD;
    
    uLL     = (LL - y)  ./ ystd;
    uDL     = (LoD - y) ./ ystd;
        
    % Call normal fsres for the ones above the detection limit
    sres = fsres( sy, yexp, ystd, fiterrors_correction_factor );

    % Overwrite selected sensitivities
    % res was: sqrt( - 2 * log( eps(1) + Phi( ( LoDidx - yidx ) ./ ystdidx ) - Phi( ( LLidx - yidx ) ./ ystdidx ) ) ); therefore
    %  sres = - ( phi(UDL)*diff(UDL, p) - phi(ULL)*diff(ULL, p) ) / ( res * ( Phi(UDL) - Phi(ULL) ) )
    if ( isinf( LL ) )
        % No lower concentration bound exists (for example when 
        denominator = repmat( res .* ( Phi(uDL) + eps(1) ), 1, 1, size( sy, 3 ) );
        phi_UDL = repmat( phi(uDL), 1, 1, size( sy, 3 ) );
        sres_bdl = -( phi_UDL .* diffUDL_p( y, sy, systd, ystd, LoD ) ) ./ denominator;        
    else
        denominator = repmat( res .* ( Phi(uDL) - Phi(uLL) + eps(1) ), 1, 1, size( sy, 3 ) );
        phi_UDL = repmat( phi(uDL), 1, 1, size( sy, 3 ) );
        phi_ULL = repmat( phi(uLL), 1, 1, size( sy, 3 ) );
        sres_bdl = -( phi_UDL .* diffUDL_p( y, sy, systd, ystd, LoD ) - phi_ULL .* diffULL_p( y, sy, systd, ystd, LL ) ) ./ denominator;        
    end
    
    % Overwrite the ones involving the values below the detection limit
    sres( repmat( idx, 1, 1, size( sres, 3 ) ) ) = sres_bdl( repmat( idx, 1, 1, size( sres, 3 ) ) );
end

% Note that only values above LoD get a sigma term here, since the values
% below the likelihood already have this term taken into account directly.
function [reserr,chi2err] = fres_error_LoD(yexp, ystd, add_c, LoD, LL)
    idx = yexp <= repmat( LoD, size(yexp,1), 1 );
    reserr = fres_error(ystd, add_c);
    
    % Sigma is already taken into account in the LoD error part
    reserr( idx ) = 0;
    
    chi2err = sum(reserr.^2,1) - add_c*sum(abs(reserr)>0,1);
end

function [sres,sreserr] = fsres_error_LoD(yexp, y, sy, ystd, systd, res, reserr, add_c, fiterrors_correction_factor, sres, LoD, LL)
    idx = yexp <= repmat( LoD, size(yexp,1), 1 );
    idxs = repmat( idx, 1, 1, size(sres, 3) );
       
    [sres_tmp, sreserr] = fsres_error(yexp, res, reserr, sres, ystd, systd);
    
    % The residuals below the detection limit don't need additional terms
    % and are already complete. The ones above, need to have the additional
    % terms included.
    sres(~idxs) = sres_tmp(~idxs);
    
    % The LoD ones don't have this term   
    sreserr(idxs) = zeros(size(sreserr(idxs)));
end


%% Standard least squares functions. These are used > limit of detection
function [res,chi2] = fres(yexp, y, ystd, fiterrors_correction_factor)
    res = (yexp-y)./ystd * sqrt(fiterrors_correction_factor);
    res(isnan(res)) = 0;
    chi2 = sum(res.^2,1);
end

function sres = fsres(sy, yexp, ystd, fiterrors_correction_factor)
sres = NaN(size(sy));
    for ip=1:size(sy,3)
        sres(:,:,ip) = - sy(:,:,ip) ./ ystd * sqrt(fiterrors_correction_factor);
        tmp = sres(:,:,ip);
        tmp(isnan(yexp)) = 0;
        tmp(isinf(yexp)) = 0;
        sres(:,:,ip) = tmp;
    end
end

function [reserr,chi2err] = fres_error(ystd, add_c)

    reserr = 2.0*log(ystd) + add_c;    
    reserr(isnan(ystd)) = 0;
    
    if(sum(reserr(:) < 0)>0)
        error('arCalcRes/fres_error: error residual too small.');
    else 
        reserr = sqrt(reserr);
        chi2err = sum(reserr.^2,1) - add_c*sum(abs(reserr)>0,1);  
    end    
end

function [sres,sreserr] = fsres_error(yexp, res, reserr, sres, ystd, systd)
    sreserr = NaN(size(sres));
    for ip=1:size(sres,3)
        sres(:,:,ip) = sres(:,:,ip) - (systd(:,:,ip) .* res ./ ystd);
        sreserr(:,:,ip) = systd(:,:,ip) ./ (reserr.*ystd);

        tmp = sres(:,:,ip);
        tmp(isnan(yexp)) = 0;
        sres(:,:,ip) = tmp;

        tmp = sreserr(:,:,ip);
        tmp(isnan(yexp)) = 0;
        sreserr(:,:,ip) = tmp;
    end
end