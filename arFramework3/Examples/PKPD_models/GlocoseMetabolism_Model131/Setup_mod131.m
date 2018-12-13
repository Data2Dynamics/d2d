%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000131#Overview)
% initialize model
arInit
arLoadModel('mod131');
arLoadData('SimuInsulinAction','mod131')
arCompileAll;


% initial values from DDmore (ERROR_prop is zero)
ar.qLog10(3) = 0;
ar.qFit(3) = 0;
varInitWOinit = [log10(150),log10(0.004),0,log10(0.03),log10(0.3),log10(0.03),log10(1),log10(0.3),log10(0.02),log10(0.4),log10(0.5)];
ar.p = varInitWOinit;


