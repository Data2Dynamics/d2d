%  [chi2_out, xSim, it, xSens_tmp, xSim2_out] = arPPL_GetChi2(t_tmp, doPPL_stuff, general_struct, stepsize, tmp_xFit, varargin)
% 
% used by 
%   arIntegratePredBand.m
%   arPPL.m

function [chi2_out, xSim, it, xSens_tmp, xSim2_out] = arPPL_GetChi2(t_tmp, doPPL_stuff, general_struct, stepsize, tmp_xFit, varargin)

global ar;

%fill temporary variables
pReset = ar.p;
data_cond = general_struct.data_cond;
x_y = general_struct.x_y;
m=general_struct.m;
c=general_struct.c; 
jx=general_struct.jx; 
takeY=general_struct.takeY; 
qLog10=ar.ppl.qLog10;
p_chi2 = ar.p(ar.qFit==1);
xstd = ar.ppl.options.xstd;

if(~isempty(varargin)>0)
    if(length(varargin{1})>1)
        p_chi2 = varargin{1};   
    end
end
if(size(p_chi2,1)~=1)
    p_chi2 = p_chi2';
end

if(doPPL_stuff)
    arLink(true,t_tmp+stepsize,takeY,jx, c, m,tmp_xFit,xstd);    
    arCalcMerit(ar.config.useSensis, p_chi2,1);
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp-stepsize));
    if(takeY)
        if(length(find(ar.model(m).(data_cond)(c).tExp==t_tmp+stepsize))>1)
            it = it+1;
        end
        xSens_tmp = -ar.sres(ar.ppl.resi_tmp,ar.qFit==1) * xstd;  
        xSim2 = ar.model(m).(data_cond)(c).yExpSimu(it,jx);
    else
        sxSim = zeros(1,length(ar.p));
        sx_tmp = arTrafoParameters(ar.model(m).(data_cond)(c).sxExpSimu,m,c,general_struct.takeY);
        sxSim(ar.model(m).(data_cond)(c).pLink) = ...
            squeeze(sx_tmp(it,jx,:))';
%         for j10=find(ar.qLog10==1)
%             sxSim(j10) = sxSim(j10) * 10.^ar.p(j10) * log(10);
%         end
        
        if(qLog10)
            sxSim = sxSim / 10^xSim / log(10);
        end
        xSens_tmp = sxSim(ar.qFit==1);
        xSim2 = ar.model(m).(data_cond)(c).xExpSimu(it,jx);
    end
end

arLink(true,t_tmp,takeY,jx, c, m,tmp_xFit,xstd);

try
    arCalcMerit(ar.config.useSensis, p_chi2,1)
catch
    fprintf('arCalcMerit doesnt work, exiting!\n');
    arLink(true,ar.model(m).(data_cond)(c).tExp(1),takeY,jx, c, m,NaN);
    ar.p = pReset;
    chi2_out=NaN;
    return;
end
[~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));

if(takeY && length(find(ar.model(m).(data_cond)(c).tExp==t_tmp))>1)
    it = it+1;
end
xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);

if(qLog10)
    xSim = log10(xSim);
end

if(doPPL_stuff)
    xSim2_out = (xSim2 - xSim)/abs(stepsize);
end

if(~takeY)
    chi2_out = arGetMerit('chi2') + ((tmp_xFit-xSim)/xstd).^2;
else
    chi2_out = arGetMerit('chi2');
end

