%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   Setup                                   %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Setup file for Final PanRTK model. Loads model file, all data and
%  compiles it. Then, some fixed parameters are set and the best parameter
%  estimate is loaded. 
%  Afterwards, either plots can be created, or validation fits for random
%  receptor/surface levels or model response for ligand stimulation wrt
%  basal level for a specific cell line

arInit;
addpath('../PanRTK_commonScripts')
%Load model and data files
arLoadModel('RTKmodel_Met_qFACS')
arLoadData('RTKdata_332M_qFACS',1,'csv');
arLoadData('RTKdata_BxPc3_qFACS',1,'csv');
arLoadData('RTKdata_A431_qFACS',1,'csv');
arLoadData('RTKdata_BT20_qFACS',1,'csv');
arLoadData('RTKdata_ADRr_qFACS',1,'csv');
arLoadData('RTKdata_IGROV1_qFACS',1,'csv');
arLoadData('RTKdata_ACHN_Met',1,'csv');
arLoadData('RTKdata_ADRr_BTC',1,'csv');

arCompileAll

%set some fixed parameters and modeling options
arSetParsPattern('basal_', -1, 1, 1, -5, 6);
arSetParsPattern('scale_', -1, 1, 1, -4, 6);
arSetParsPattern('sd_', -1, 2, 1, -3, 6);
arSetParsPattern('relto_', -1, 1, 1, -2.5, 2.5);
arSetPars('scale_Ligand',2,1,1,1,6);
arSetParsPattern('offset_',-1,1,1,-5,5);
arSetParsPattern('offset_pERK',-7,2,1,-8,-6);
arSetParsPattern('offset_pS6K1',-7,2,1,-8,-6);
arSetParsPattern('offset_pS6',-7,2,1,-8,-6);

ar.config.fiterrors=0;
ar.config.maxsteps=1.e4;

arSetPars({'init_EGFR_EGF','init_ErbB3_HRG','init_IGF1R_IGF1','init_Met_HGF','init_EGFR_BTC'},zeros(1,5),ones(1,5)*2,zeros(1,5),ones(1,5)*-5,ones(1,5)*3);
arSetParsPattern('sd_FACS_',-1.2,2,1,-5,3);
arSetPars({'EGF_kD','HRG_kD','IGF1_kD','HGF_kD'},[log10(1) log10(0.05) log10(0.3) log10(0.3)],ones(1,4)*2,ones(1,4),ones(1,4)*-5,ones(1,4)*3);
%Load pars of best fit and do the plotting
arLoadPars('Final_Model')
arSave('Model_Final')
arCalcMerit
