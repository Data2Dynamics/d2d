% addpath('../TransientFunction_library/') % functions for fitting the transient model (e.g. fitting both signums)

nLHS = 500;

%% LDH
arInit
arLoadModel('TransientDoseFunction')
% arLoadData('LDH',1)
arLoadData('LDH_WT',1)
arLoadData('LDH_KO',1)
arCompileAll(true)

% Initialize_FitTransient % setting bounds, identification of parameter-types from name

% arLoadPars(1)

arSetPars('KD_toffset_TF_sus',log10(2), 1,1,-5,3);
arSetPars('Max_toffset_TF_sus',10, 1,0,-1,100); % unlog
arSetPars('hill_toffset_TF_sus',1, 1,0,1,5);
arSetPars('hill_timescale_sus',1, 1,0,1,5);
arSetPars('hill_amp_sus',1, 1,0,1,5);

arFitLHS(nLHS) % multi-start fitting
% arFitTransient % single fit
arQplot('y')
arPlot
arPrint

%% 
arSave('LDH')

arPLEInit
ple
pleExtend
pleSmooth
arSave('current')

%%
% arLoad('20220815T133125_PLEs_Gut')
PlotHillCurves(2)
print -dpng HillCurves_LDH

PlotHillCurves(3)
print -dpng HillCurves2_LDH

PlotHillCurves_3D
PrintAllToPng(1:2,'HillCurves_LDH_3D')

% PLE 
PlotProfile('fold_as','max. amplitude ratio KO/WT')
title('LDH')
print -dpng Profile_fold_Amax_LDH

PlotProfile('fold_ts','velocity ratio KO/WT')
title('LDH')
print -dpng Profile_fold_timescale_LDH

PlotProfile('fold_to','delay ratio KO/WT')
title('LDH')
print -dpng Profile_fold_toffset_LDH

PlotFits
suptitle('LDH')
print -dpng PlotFits_LDH


Plot3D
print -dpng 3D_LDH
% arSave('LDH_foldAufHill_LHS500_gut')

%% IL1b
arInit
arLoadModel('TransientDoseFunction')
% arLoadData('LDH',1)
arLoadData('IL1b_WT',1)
arLoadData('IL1b_KO',1)
arCompileAll(true)

% Initialize_FitTransient % setting bounds, identification of parameter-types from name

% arLoadPars(1)

arSetPars('KD_toffset_TF_sus',log10(2), 1,1,-5,3);
arSetPars('Max_toffset_TF_sus',10, 1,0,-1,100); % unlog
arSetPars('hill_toffset_TF_sus',1, 1,0,1,5);
arSetPars('hill_timescale_sus',1, 1,0,1,5);
arSetPars('hill_amp_sus',1, 1,0,1,5);

arFitLHS(nLHS) % multi-start fitting
% arFitTransient % single fit
arQplot('y')
arPlot
arPrint

%%
arSave('IL1b')

arPLEInit
ple
pleExtend
arSave('current')


%% without t_offset
% indFix = [arPrint('fold_to'),arPrint('_toffset_')]
% ar.qFit(indFix) = 0;

arFit
arSave('IL1b_toffsetFixed')

arPLEInit
ple
pleExtend
pleSmooth
arSave('current')

%% Fits IL1b
% arLoad('20220815T150340_IL1b_toffsetFixed')

PlotFits
suptitle('IL1b')
print -dpng PlotFits_IL1b

PlotProfile('fold_as','max. amplitude ratio KO/WT')
title('IL1b')
print -dpng Profile_fold_Amax_IL1b

PlotProfile('fold_ts','velocity ratio KO/WT')
title('IL1b')
print -dpng Profile_fold_timescale_IL1b

PlotProfile('fold_to','delay ratio KO/WT')
title('IL1b')
print -dpng Profile_fold_toffset_IL1b


PlotHillCurves(2)
print -dpng HillCurves_IL1b

PlotHillCurves(3)
print -dpng HillCurves2_IL1b

PlotHillCurves_3D
PrintAllToPng(1:2,'HillCurves_IL1b_3D')

PlotHillCurves_3D('1./((Max_timescale_sus*KD_timescale_sus^(hill_timescale_sus))./(KD_timescale_sus^(hill_timescale_sus)+dosesFine.^(hill_timescale_sus))*(1-isKO) + isKO*fold_ts*(Max_timescale_sus*KD_timescale_sus^(hill_timescale_sus))./(KD_timescale_sus^(hill_timescale_sus)+dosesFine.^(hill_timescale_sus)));',...
    'Velocity');
PlotHillCurves_3D('((Max_toffset_TF_sus*KD_toffset_TF_sus^(hill_toffset_TF_sus))./(KD_toffset_TF_sus^(hill_toffset_TF_sus)+dosesFine.^(hill_toffset_TF_sus))*(1-isKO) + isKO*fold_to*(Max_toffset_TF_sus*KD_toffset_TF_sus^(hill_toffset_TF_sus))./(KD_toffset_TF_sus^(hill_toffset_TF_sus)+dosesFine.^(hill_toffset_TF_sus)));',...
    'Max. Delay');

Plot3D
print -dpng 3D_IL1b

% arSave('IL1b_foldAufHill_LHS500_gut')
