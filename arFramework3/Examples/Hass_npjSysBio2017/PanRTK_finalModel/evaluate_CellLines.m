%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   evaluate_CellLines                      %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Go through all calibration/validation cell lines to collect their
%  goodness of fit and do n calculations with random assigned
%  receptor/surface levels
%
%
function evaluate_CellLines(nRnd, doMore)
global ar
if(~exist('doMore','var') || isempty(doMore))
    doMore = 0;
end
%validation cell lines
ucelllines = {'332M', 'BxPc3','A431','BT20','ADRr','IGROV1','BT474','ACHN' 'MDAMB231'};%,
if(doMore && nRnd+length(ucelllines) <= length(ar.val.chi2_final))
   error('caliculated already that amount of random qFACS'); 
end
nr_val = length(ucelllines)-6;
if(~doMore)
    if(length(ar.model(1).plot)>6)
       ar.model(1).plot(7:end)=[]; 
    end
    ar.val.nr_val = nr_val;
    ar.val.plots=length(ar.model(1).plot);
end
ar.val.nRnd = nRnd;
ar.val.chi2_detailed = NaN(1,length(ar.model(1).yNames)*10);
ar.val.ndata_detailed = NaN(1,length(ar.model(1).yNames)*10);
ar.val.chi2_qFit = NaN(1,length(ar.model(1).yNames)*10);
starting = 1;
if(doMore)
    starting = length(ar.val.chi2_final)+1;
end
for i=starting:(ar.val.plots + ar.val.nr_val + nRnd)
    ndata_tmp = 0;
    chi2_tmp = 0;
    %cell lines used for calibration
    if(i<=ar.val.plots)% && sum(find(~cellfun(@isempty,strfind(ucelllines,ar.model(1).plot(i).name))))==0)
        for jd=ar.model(1).plot(i).dLink
            ndata_tmp = ndata_tmp + sum(ar.model(1).data(jd).ndata(ar.model(1).data(jd).qFit==1));
            chi2_tmp = chi2_tmp + sum(ar.model(1).data(jd).chi2(ar.model(1).data(jd).qFit==1));
        end
    %twice looping through validation cell lines
    else
        load_cell = sprintf('RTKdata_%s_qFACS',ucelllines{mod(i-1,length(ucelllines))+1});
        LoadABData(load_cell,'csv');
        %in second loop, random qFACS; yExp is log10, qFACS is normalized
        %by 10000, thus 0 corresponds to qFACS==1.e4, 2 is 1.e6
        if(i>ar.val.plots+ar.val.nr_val)
            which_cell=randi([1 size(ar.nr_cells,2)],1);
            ar.model(1).data(1).yExp(1,1:5) = comp_FACS(ar.nr_cells([1:4 6],which_cell)');
            %get Random qFACS instead of the ones from cell lines?
            %ar.model(1).data(1).yExp(1,1:5) = [rand*2.5 rand*2.5 (rand-0.3)*2.2 (rand-0.3)*2.2 (rand-0.3)*2.2];            
        end
        fitLoaded
        for jd=ar.model(1).plot(end).dLink
            ndata_tmp = ndata_tmp + sum(ar.model(1).data(jd).ndata(ar.model(1).data(jd).qFit==1));
            chi2_tmp = chi2_tmp + sum(ar.model(1).data(jd).chi2(ar.model(1).data(jd).qFit==1));
            
            ar.val.chi2_detailed(i,(jd-1)*length(ar.model(1).yNames)+1:jd*length(ar.model(1).yNames)) = ar.model(1).data(jd).chi2;
            ar.val.chi2_qFit(i,(jd-1)*length(ar.model(1).yNames)+1:jd*length(ar.model(1).yNames)) = ar.model(1).data(jd).qFit==1;
            ar.val.ndata_detailed(i,(jd-1)*length(ar.model(1).yNames)+1:jd*length(ar.model(1).yNames)) = ar.model(1).data(jd).ndata;

        end       
        dose_plots(1,length(ar.model(1).plot));
        
    end
    
    ar.val.chi2(i) = chi2_tmp;
    ar.val.ndata(i) = ndata_tmp;
    ar.model(1).plot(i).chi2 = chi2_tmp;
    ar.model(1).plot(i).ndata = ndata_tmp;
    ar.val.chi2_final(i) = chi2_tmp/ndata_tmp;
    
    
end

ar.val.mean_cal = nanmean(ar.val.chi2_final(1:ar.val.plots));
ar.val.mean_val = nanmean(ar.val.chi2_final(ar.val.plots+1:ar.val.plots+ar.val.nr_val));
ar.val.mean_val_rnd = nanmean(ar.val.chi2_final(ar.val.plots+ar.val.nr_val+1:end));
ar.val.std_cal = nanstd(ar.val.chi2_final(1:ar.val.plots));
ar.val.std_val = nanstd(ar.val.chi2_final(ar.val.plots+1:ar.val.plots+ar.val.nr_val));
ar.val.std_val_rnd = nanstd(ar.val.chi2_final(ar.val.plots+ar.val.nr_val+1:end));


function qFACS = comp_FACS(RNA)
    qFACS(1,1) = log10(sinh(asinh(RNA(1))*0.964 + 8.9229)/10000);
    qFACS(1,2) = log10(sinh(asinh(RNA(2))*0.8282 + 8.0359)/10000);
    qFACS(1,3) = log10(sinh(asinh(RNA(3))*0.3445 + 9.1685)/10000);
    qFACS(1,4) = log10(sinh(asinh(RNA(4))*0.4757 + 9.47)/10000);
    qFACS(1,5) = log10(sinh(asinh(RNA(5))*0.6418 + 8.3577)/10000);