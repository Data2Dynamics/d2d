% Plot sensitivities of models and datasets
%
% arPlotS

function arPlotS

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

arSimu(true, true);
arSimu(true, false);

arPlotSY;
arPlotSYSTD;
arPlotSX;
