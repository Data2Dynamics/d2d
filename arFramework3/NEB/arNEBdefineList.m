function arNEBDefineList
% makes list with #define expressions for _NEB_XXX.def template

global ar

plabelsfit = ar.pLabel(ar.qFit==1);

for i = 1:length(plabelsfit)

fprintf(['#define ' plabelsfit{i} ' ' plabelsfit{i} '_NEB_XXX \n'])

end
