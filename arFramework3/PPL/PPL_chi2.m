%  [chi2_out, xSim, xSim2, xSim3, it] = PPL_chi2(t_tmp, doPPL_stuff, m, c, jx, takeY, qLog10, doPPL, stepsize, tmp_xFit, xstd, varargin)
% 
% used by 
%   doPPL.m
%   PPL_int.m
%   VPL_int.m

function [chi2_out, xSim, xSim2, xSim3, it] = PPL_chi2(t_tmp, doPPL_stuff, m, c, jx, takeY, qLog10, doPPL, stepsize, tmp_xFit, xstd, varargin)

global ar;
xSim2=NaN;
xSim3=NaN;
RHS_t = NaN;
pReset = ar.p;
if(~exist('doPPL_stuff','var'))
    doPPL_stuff = false;
end
if(takeY)
    data_cond = 'data';
    x_y = 'y';
else
    data_cond = 'condition';
    x_y = 'x';
end
p_chi2 = ar.p(ar.qFit==1);
if(~isempty(varargin)>0)
    if(length(varargin{1})>1)
        p_chi2 = varargin{1};
    else
        RHS_t =  varargin{1};
    end
end
if(size(p_chi2,1)~=1)
    p_chi2 = p_chi2';
end
if(nargout>2)
    get_sensi=true;
else
    get_sensi=true;
end

if(doPPL_stuff)
    arLink(true,t_tmp+stepsize,takeY,jx, c, m,tmp_xFit,xstd);    
    arChi2(ar.config.useSensis, p_chi2,1);
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp-stepsize));
    if(takeY)
        if(length(find(ar.model(m).(data_cond)(c).tExp==t_tmp+stepsize))>1)
            it = it+1;
        end
        ar.model(m).(data_cond)(c).ppl.xSens_tmp = -ar.sres(ar.ppl.resi_tmp,ar.qFit==1) * xstd;
        xSim2 = ar.model(m).(data_cond)(c).yExpSimu(it,jx);
        
    else
        sxSim = zeros(1,length(ar.p));
        sxSim(ar.model(m).(data_cond)(c).pLink) = ...
            squeeze(ar.model(m).(data_cond)(c).sxExpSimu(it,jx,:))';
        for j10=find(ar.qLog10==1)
            sxSim(j10) = sxSim(j10) * 10.^ar.p(j10) * log(10);
        end
        
        if(qLog10)
            sxSim = sxSim / 10^xSim / log(10);
        end
        ar.model(m).(data_cond)(c).ppl.xSens_tmp = sxSim(ar.qFit==1);
        xSim2 = ar.model(m).(data_cond)(c).xExpSimu(it,jx);
    end
    arLink(true,t_tmp-1.e-3,takeY,jx, c, m,tmp_xFit,xstd);
    arChi2(false, p_chi2,1);
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp+1.e-3));
    xSim3 = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);

    if(~isnan(RHS_t))
        arLink(true,t_tmp+RHS_t,takeY,jx, c, m,tmp_xFit,xstd);
        arChi2(false, p_chi2,1);
        [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp-RHS_t));
        ar.ppl.fRHS_ppl = ar.model(m).(data_cond)(c).dxdts(it,jx);
    end
end
arLink(true,t_tmp,takeY,jx, c, m,tmp_xFit,xstd);

try
    arChi2(get_sensi, p_chi2,1)
catch
    fprintf('arChi2 doesnt work, exiting for direction %i!\n',dir);
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
    %x_orig=log10(x_orig);
end

if(doPPL_stuff)
    ar.ppl.xSim2 = (xSim2 - xSim)/abs(stepsize);
end

res = ar.res;
if(~takeY && ~doPPL)
    res(end+1) = (tmp_xFit-ar.model(m).condition(c).xExpSimu(it,jx))/xstd;
elseif(takeY && doPPL)
    res(ar.ppl.resi_tmp) = [];
end
chi2_out = nansum(res.^2);

