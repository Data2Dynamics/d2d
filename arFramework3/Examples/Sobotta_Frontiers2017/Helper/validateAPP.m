validation = arFindData('validation', 'verbose');

% Lock all parameters
ar.qFit = 0 * ar.qFit;

% Unlock the validation parameters
ar.qFit(arFindPar('validation')) = 1;

arDisableData( 'all', 1 );
arEnableData( arFindData('validation', 'names'), 1 );