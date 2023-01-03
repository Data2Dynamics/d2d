
function TransPars = arTransientPars

global ar

y = ar.model(1).data.yFineSimu;
t = ar.model(1).data.tFine;
yname = ar.model.data.y;

ar.config.fiterrors=1;  % ignore errors of biomodel data file

for i = 1:size(y,2)

    dat.tExp = t(:,1);
    if i>1 && size(t,2)>1
        dat.tExp = t(:,i);
    end
    dat.yExp = y(:,i);
    dat.ystd = nan(size(dat.yExp));

    res = arFitTransientFunction2(dat,[pwd filesep 'RealisticDesign' filesep 'TransientFit_' yname{i}]);
    %arPlotY;
    if i==1
        TransPars = table('Size',[size(y,2) length(ar.p)],'VariableTypes',repmat("double",1,length(ar.p)),'VariableNames',ar.pLabel);
    end
    TransPars{i,:} = res.pRescaled;
end
% system('ps2pdf RealisticDesign/TransientFit.ps');
% delete 'RealisticDesign/TransientFit.ps'
% taus = TransPars(:,contains(TransPars.Properties.VariableNames,'time') | contains(TransPars.Properties.VariableNames,'toffset'));

writetable([array2table(yname'), TransPars],'RealisticDesign/TransPars.txt');
fprintf('Time parameters assessed by fitting transient function.\n');