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
if(takeY)
    
    if(doPPL_stuff)
        arLink(true,t_tmp+stepsize,true,jx, c, m,tmp_xFit,xstd);
        arChi2(ar.config.useSensis, p_chi2);
        [~,it] = min(abs(ar.model(m).data(c).tExp-t_tmp-stepsize));
        if(length(find(ar.model(m).data(c).tExp==t_tmp+stepsize))>1)
            it = it+1;
        end
        ar.model(m).data(c).ppl.xSens_tmp = -ar.sres(ar.ppl.resi_tmp,ar.qFit==1) * xstd;
        %ar.sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,ar.qFit==1)
        xSim2 = ar.model(m).data(c).yExpSimu(it,jx);
        
        arLink(true,t_tmp-1.e-3,true,jx, c, m,tmp_xFit,xstd);
        arChi2(false, p_chi2);
        [~,it] = min(abs(ar.model(m).data(c).tExp-t_tmp+1.e-3));
        xSim3 = ar.model(m).data(c).yExpSimu(it,jx);
        
        if(~isnan(RHS_t))
            arLink(true,t_tmp+RHS_t,true,jx, c, m,tmp_xFit,xstd);
            arChi2(false, p_chi2);
            [~,it] = min(abs(ar.model(m).data(c).tExp-t_tmp-RHS_t));
            ar.ppl.fRHS_ppl = ar.model(m).data(c).dxdts(it,jx);
        end
    end
    
    arLink(true,t_tmp,true,jx, c, m,tmp_xFit,xstd);
    try
        arChi2(get_sensi, p_chi2)
    catch
        fprintf('arChi2 doesnt work, exiting for direction %i!\n',dir);
        arLink(true,ar.model(m).data(c).tExp(1),true,jx, c, m,NaN);
        ar.p = pReset;
        chi2_out=NaN;
        return;
    end
    [~,it] = min(abs(ar.model(m).data(c).tExp-t_tmp));
    
    if(length(find(ar.model(m).data(c).tExp==t_tmp))>1)
        it = it+1;
    end
    xSim = ar.model(m).data(c).yExpSimu(it,jx);
    
else
    if(doPPL_stuff)
        arLink(true,t_tmp+stepsize);
        arChi2(ar.config.useSensis, p_chi2);
        [~,it] = min(abs(ar.model(m).condition(c).tExp-t_tmp-stepsize));
        sxSim = zeros(1,length(ar.p));
        sxSim(ar.model(m).condition(c).pLink) = ...
            squeeze(ar.model(m).condition(c).sxExpSimu(it,jx,:))';
        for j10=find(ar.qLog10==1)
            sxSim(j10) = sxSim(j10) * 10.^ar.p(j10) * log(10);
        end
        
        if(qLog10)
            sxSim = sxSim / 10^xSim / log(10);
        end
        ar.model(m).condition(c).ppl.xSens_tmp = sxSim(ar.qFit==1);
        xSim2 = ar.model(m).condition(c).xExpSimu(it,jx);
        
        arLink(true,t_tmp-1.e-3);
        arChi2(false, p_chi2);
        [~,it] = min(abs(ar.model(m).condition(c).tExp-t_tmp+1.e-3));
        xSim3 = ar.model(m).condition(c).xExpSimu(it,jx);
        
        if(~isnan(RHS_t))
            arLink(true,t_tmp+RHS_t);
            arChi2(false, p_chi2);
            [~,it] = min(abs(ar.model(m).condition(c).tExp-t_tmp-RHS_t));
            ar.ppl.fRHS_ppl = ar.model(m).condition(c).dxdts(it,jx);
        end
        
    end
    arLink(true,t_tmp);
    try
        arChi2(get_sensi, p_chi2)
    catch
        fprintf('arChi2 doesnt work, exiting for t %f !\n',t_tmp);
        ar.p = pReset;
        chi2_out=NaN;
        return;
    end
    [~,it] = min(abs(ar.model(m).condition(c).tExp-t_tmp));
    xSim = ar.model(m).condition(c).xExpSimu(it,jx);
    
end
if(qLog10)
    xSim = log10(xSim);
    %x_orig=log10(x_orig);
end

if(doPPL_stuff)
    ar.ppl.xSim2 = (xSim2 - xSim)/abs(stepsize);
end
%elseif(doPPL_stuff)
%ar.ppl.xSim2 = (xSim - xSim3)/1.e-4;
%end
%res = [ar.res ar.constr];
res = ar.res;
if(~takeY && ~doPPL)
    res(end+1) = (tmp_xFit-ar.model(m).condition(c).xExpSimu(it,jx))/xstd;
elseif(takeY && doPPL)
    res(ar.ppl.resi_tmp) = [];
end
chi2_out = nansum(res.^2);

