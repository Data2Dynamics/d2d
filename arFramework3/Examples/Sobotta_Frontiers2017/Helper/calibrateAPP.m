dataSets;

% Fit core along with APP
incremental = 0;

% Lock all pars
if (~incremental)
    ar.qFit = zeros(size(ar.qFit));
end
    
% Unlock APP calibration parameters
ar.qFit(calibrationPars) = 1;

% Disable all data
if (~incremental)
    ar = arDisableData( 'all' );
end

% Enable only the APP calibration data
ar = arEnableData( app_calibrationdataSets );