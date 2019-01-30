% [dps, dx_out, gamma_tmp] = arIntStepPPL(t_tmp, general_struct, stepsize, gamma_tmp)
% 
% t_tmp
% ppl_general_struct    PPL result struct
%                       model index, condition index etc are taken from
%                       this struct
% stepsize
% gamma_tmp
% 
% dps
% dx_out
% gamma_tmp
% 
% This function calls arPPL_GetChi2
% 
% Written by Helge, tried to be documented by Clemens.

function [dps, dx_out, gamma_tmp] = arIntStepPPL(t_tmp, general_struct, stepsize, gamma_tmp)
global ar;

m=general_struct.m;
c=general_struct.c;
jx=general_struct.jx; 
takeY=general_struct.takeY;
qLog10 = ar.ppl.qLog10;
xstd = ar.ppl.options.xstd;
%only needed if RK4 is chosen
xFit_tmp = ar.ppl.xFit_tmp;
p_tmp = general_struct.pReset;

[grad, H, dy_pred, xSim, grad_y, fRHS] = get_stuff(t_tmp,stepsize);
dx_0 = get_int();
region_fac_p = get_fac(dx_0(2:end));       
dx_0(2:end) = dx_0(2:end).*region_fac_p';

%Just do Euler integration
%if((min(region_fac_p)==1 || ~ar.ppl.options.fineInt)
dx_out = dx_0;
dps = dx_0(2:end);
return;

%Uncomment if statement and comment out return to obtain RK4 integration
%else
   ar.p(ar.qFit==1) = ar.p(ar.qFit==1) + dx_0(2:end)'*stepsize/2;
   ar.ppl.xFit_tmp = ar.ppl.xFit_tmp + dx_0(1)*stepsize/2;
%end
[grad, H, dy_pred, xSim, grad_y, fRHS] = get_stuff(t_tmp+stepsize/2,stepsize);
dx_A = get_int();
region_fac_p = get_fac(dx_A(2:end));
if(max(region_fac_p) < 0.1)
    region_fac_p = zeros(1,length(dx_0(2:end)));
end 
dx_A(2:end) = dx_A(2:end).*region_fac_p';
ar.p(ar.qFit==1)=p_tmp(ar.qFit==1) + dx_A(2:end)'*stepsize/2;
ar.ppl.xFit_tmp = xFit_tmp + dx_A(1)*stepsize/2;
ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);        
[grad, H, dy_pred, xSim, grad_y, fRHS] = get_stuff(t_tmp+stepsize/2,stepsize);
dx_B = get_int();
region_fac_p = get_fac(dx_B(2:end));
if(max(region_fac_p) < 0.1)
    region_fac_p = zeros(1,length(dx_0(2:end)));
end 
dx_B(2:end) = dx_B(2:end).*region_fac_p';
ar.p(ar.qFit==1)=p_tmp(ar.qFit==1) + dx_B(2:end)'*stepsize;
ar.ppl.xFit_tmp = xFit_tmp + dx_B(1)*stepsize;
ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);        
[grad, H, dy_pred, xSim, grad_y, fRHS] = get_stuff(t_tmp+stepsize,stepsize);
dx_C = get_int();
ar.p(ar.qFit==1)=p_tmp(ar.qFit==1);
ar.ppl.xFit_tmp = xFit_tmp;
region_fac_p = get_fac(dx_C(2:end));
if(max(region_fac_p) < 0.1)
    region_fac_p = zeros(1,length(dx_0(2:end)));
end 
dx_C(2:end) = dx_C(2:end).*region_fac_p';
dx_out = 1/6*(dx_0+2*(dx_A+dx_B)+dx_C);
dps = dx_out(2:end);


    function [grad, H, dy_pred, xSim, grad_y, fRHS] = get_stuff(t_get, step)
        %The arPPL_GetChi2 function computes the required sensitivities and
        %values of RHS. resi_tmp is set in arLink
        [~, xSim, it, xSens_tmp, fRHS] = arPPL_GetChi2(t_get, true, general_struct, step, ar.ppl.xFit_tmp);
        res = ar.res;
        if(~takeY)
           res(end+1) = (ar.ppl.xFit_tmp - xSim).^2 ./ xstd.^2; 
        end
        if(ar.config.useSensis)
            sres = ar.sres(:, ar.qFit==1);               
        %                 if(~isempty(ar.sconstr))
        %                     sres = [ar.sres(:, ar.qFit==1); ar.sconstr(:, ar.qFit==1)];
        %                 end
            if(~takeY)        
                sxSim = zeros(1,length(ar.p));
                sx_tmp = arTrafoParameters(ar.model(m).condition(c).sxExpSimu,m,c,false);
                sxSim(ar.model(m).condition(c).pLink) = ...
                    squeeze(sx_tmp(it,jx,:))';
%                 for j=find(ar.qLog10==1)
%                     sxSim(j) = sxSim(j) * 10.^ar.p(j) * log(10);
%                 end

                if(qLog10)
                    sxSim = sxSim / 10^xSim / log(10);
                end

                sres(end+1,:) = -sxSim(ar.qFit==1) / xstd;
            end
        end       

        grad = 2*res*sres;
        H = 2*(sres')*sres;
        dy_pred = zeros(1,length(ar.p));
        if(takeY)            
            dy_pred(ar.qFit==1) = (xSens_tmp + sres(ar.ppl.resi_tmp,:) * xstd )/(abs(step) * xstd^2); %sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:)
            grad_y = -sres(ar.ppl.resi_tmp,:);
        else
            dy_pred(ar.qFit==1) =(xSens_tmp + sres(end,:) * xstd) / (abs(step) * xstd^2); 
            grad_y = -sres(end,:);
            if(qLog10)
                fRHS = fRHS / log(10) / 10^xSim;
                dy_pred = dy_pred / 10^xSim / log(10);
            end
        end
        dy_pred = dy_pred(ar.qFit==1);
    end




    function dxs = get_int()
       Mty=zeros(length(ar.p(ar.qFit==1))+1,length(ar.p(ar.qFit==1))+1);
        Mty(1,1)=1;
        Mty(2:end,1)=0;
        Mty(1,2:end)=0;
        Mty(2:end,2:end)=H;% + diag(ones(length(ar.p(ar.qFit==1)),1))*1.e-4;   
        if((takeY && ar.model(m).data(c).logfitting(jx)==1) || qLog10)
            Mty(2:end, 2:end) = Mty(2:end, 2:end) + (ar.ppl.xFit_tmp - xSim)*log(10)*H;
        end
        %magic formula                
        if((takeY && ar.model(m).data(c).logfitting(jx)==1) || qLog10)
            dx_woGamma = pinv(Mty,1.e-6)*[fRHS;(ar.ppl.xFit_tmp-xSim)*2.*dy_pred' - 2*(ar.ppl.xFit_tmp - xSim)./xstd./(10.^xSim)*grad_y'*fRHS];            
        else
            dx_woGamma = pinv(Mty,1.e-6)*[fRHS;(ar.ppl.xFit_tmp-xSim)*2.*dy_pred'];
        end             
        dx_Gamma = pinv(Mty(2:end,2:end),1.e-6)*(- gamma_tmp.*grad');
       
        %Check if integration without gamma correction is sufficient
        region_fac_tmp = get_fac(dx_woGamma(2:end));  
        [chi2, ~] = arPPL_GetChi2(t_tmp + stepsize,false, general_struct, stepsize, ar.ppl.xFit_tmp + dx_woGamma(1)*stepsize,ar.p(ar.qFit==1)+stepsize*dx_woGamma(2:end)'*region_fac_tmp);
        if((chi2 - ar.ppl.chi2_threshold) < 0.2 && (chi2 - ar.ppl.chi2_threshold) > -0.2)
            dxs = dx_woGamma;
            return;
        end
        %Gamma control via norm of Delta p vectors
%         if(norm(dx_woGamma(2:end)) / norm(dx_Gamma) < 2)
%             gamma_tmp = gamma_tmp/2;            
%         end
%         
%         if(norm(dx_woGamma(2:end)) / norm(dx_Gamma) > 4)% && gamma_tmp < gamma_norm)
%             gamma_tmp = 2*gamma_tmp;
%             if(gamma_tmp > gamma_norm)
%                 gamma_tmp = gamma_norm;
%             end
%         end
        if((takeY && ar.model(m).data(c).logfitting(jx)==1) || qLog10)
            dxs=pinv(Mty,1.e-6)*[fRHS;((ar.ppl.xFit_tmp-xSim)*2.*dy_pred' - 2*(ar.ppl.xFit_tmp - xSim)./xstd./(10.^xSim)*grad_y'*fRHS - gamma_tmp.*grad')];            
        else
            dxs=pinv(Mty,1.e-6)*[fRHS;((ar.ppl.xFit_tmp-xSim)*2.*dy_pred' - gamma_tmp.*grad')];
        end
        
        %Adapt gamma correction factor, half gamma if calculation w/o
        %gamma was better, otherwise twice, all within 0.1 < gamma < 10
        region_fac_tmp = get_fac(dxs(2:end));                      
        [chi2_sec, ~] = arPPL_GetChi2(t_tmp + stepsize,false, general_struct, stepsize, ar.ppl.xFit_tmp + dxs(1)*stepsize,ar.p(ar.qFit==1)+stepsize*dxs(2:end)'*region_fac_tmp);
        if((chi2_sec - ar.ppl.chi2_threshold) > 0.2 || (chi2_sec - ar.ppl.chi2_threshold) < -0.2)
            if(abs(chi2-ar.ppl.chi2_threshold)<abs(chi2_sec-ar.ppl.chi2_threshold) && gamma_tmp>0.1)
                gamma_tmp = gamma_tmp/2;  
            elseif(gamma_tmp<10)
                gamma_tmp = gamma_tmp*2;  
            end
        end
    end


    function region_fac_p = get_fac(dp_tmp)
        region_fac_p=1;
            
        if(sum(ar.p(ar.qFit==1) + dp_tmp'*stepsize > ar.ub(ar.qFit==1))>0 || (sum(ar.p(ar.qFit==1) + dp_tmp'*stepsize < ar.lb(ar.qFit==1))>0))
            dps_ub = dp_tmp;
            dps_ub(dp_tmp<0) = 0;
            region_fac_p = min((ar.ub(ar.qFit==1)- ar.p(ar.qFit==1))./ (dps_ub'*stepsize));
            dps_lb = dp_tmp;
            dps_lb(dp_tmp>0) = 0;
            region_fac_p = min([region_fac_p,(ar.p(ar.qFit==1) - ar.lb(ar.qFit==1))./ (abs(dps_lb)'*stepsize)]);       
        end
        if(region_fac_p < 0.1)
            region_fac_p = 0;
        end
    end

end