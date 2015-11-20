function [dps, dx, gamma_tmp] = PPL_Int(t_tmp, m, c, jx, takeY, qLog10, stepsize, xFit_tmp, xstd, dir, gamma_tmp)
global ar;
ub_fit = ar.ub(ar.qFit==1);
lb_fit = ar.lb(ar.qFit==1);
interval=0.1;
% if(dir==1)
%     interval=0.1;
% else
%     interval=-0.1;
% end  
beta = sym('beta');
p_tmp_ppl = ar.p(ar.qFit==1);
if(size(p_tmp_ppl,1)>1)
    p_tmp_ppl = p_tmp_ppl';
end
% curr_mult = 0;
% for i=1:10
%     [fRHS, dps_A, ~] = get_ppl(t_tmp, stepsize, -stepsize/4*(i-1));
%     curr_mult = ar.model(m).condition(c).ppl.xSens_tmp * dps_A';
%     if(curr_mult>0)
%         break;
%     end
% end
% i

[fRHS, dps_A, ~] = get_ppl(t_tmp, stepsize, NaN);
dx = fRHS+0.;
dps = gamma_tmp*dps_A';

%test different dp lengths
p_ten=ar.p(ar.qFit==1)+10*dps'*stepsize;
p_ten(p_ten>ub_fit) = ub_fit(p_ten>ub_fit);
p_ten(p_ten<lb_fit) = lb_fit(p_ten<lb_fit);
[chi2_ten] = PPL_chi2(t_tmp+stepsize, 0, m, c, jx, takeY, qLog10, true, stepsize, xFit_tmp, xstd, p_ten);
if(abs(chi2_ten - ar.ppl.chi2_95 + 0.5) > 0.2)
    p_norm=ar.p(ar.qFit==1)+dps'*stepsize;
    p_norm(p_norm>ub_fit) = ub_fit(p_norm>ub_fit);
    p_norm(p_norm<lb_fit) = lb_fit(p_norm<lb_fit);
    [chi2_norm] = PPL_chi2(t_tmp+stepsize, 0, m, c, jx, takeY, qLog10, true, stepsize, xFit_tmp, xstd, p_norm);
    if(abs(chi2_norm - ar.ppl.chi2_95 + 0.5) > 0.2)
        p_tenth=ar.p(ar.qFit==1)+dps'/10*stepsize;
        p_tenth(p_tenth>ub_fit) = ub_fit(p_tenth>ub_fit);
        p_tenth(p_tenth<lb_fit) = lb_fit(p_tenth<lb_fit);
        [chi2_tenth] = PPL_chi2(t_tmp+stepsize, 0, m, c, jx, takeY, qLog10, true, stepsize, xFit_tmp, xstd, p_tenth);
        if(abs(chi2_tenth - ar.ppl.chi2_95 + 0.5) > 0.2)
            p_hunth=ar.p(ar.qFit==1)+dps'/100*stepsize;
            p_hunth(p_hunth>ub_fit) = ub_fit(p_hunth>ub_fit);
            p_hunth(p_hunth<lb_fit) = lb_fit(p_hunth<lb_fit);
            [chi2_hunth] = PPL_chi2(t_tmp+stepsize, 0, m, c, jx, takeY, qLog10, true, stepsize, xFit_tmp, xstd, p_hunth);
            if(abs(chi2_hunth - ar.ppl.chi2_95 + 0.5) > 0.2)
                dps = 0*dps;
                gamma_tmp = gamma_tmp/100;
            else
                dps = dps/100;
                gamma_tmp = gamma_tmp/100; 
            end
        else
            gamma_tmp = gamma_tmp/10;
            dps = dps/10;
        end      
    end
else
    p_hun=ar.p(ar.qFit==1)+dps'*100*stepsize;
    p_hun(p_hun>ub_fit) = ub_fit(p_hun>ub_fit);
    p_hun(p_hun<lb_fit) = lb_fit(p_hun<lb_fit);
    [chi2_hun] = PPL_chi2(t_tmp+stepsize, 0, m, c, jx, takeY, qLog10, true, stepsize, xFit_tmp, xstd, p_hun);
    if(abs(chi2_hun - ar.ppl.chi2_95 + 0.5) > 0.2)
        gamma_tmp = gamma_tmp * 10;
        dps = 10*dps;
    else
        gamma_tmp = gamma_tmp * 100;
        dps = 100*dps;
    end
end
% if(max(dps_A)==0)
%     [fRHS, dps_A] = get_ppl(t_tmp, stepsize*2, NaN);
%     dx = fRHS+0.;
%     dps = dps_A';
% end

% ar.p(ar.qFit==1) = p_tmp_ppl + dps_A;
% 
% [~, dps_B] = get_ppl(t_tmp, stepsize/2, stepsize/2);
% ar.p(ar.qFit==1) = p_tmp_ppl + dps_B;
% 
% [~, dps_C] = get_ppl(t_tmp, stepsize, stepsize/2);
% ar.p(ar.qFit==1) = p_tmp_ppl;
% dps = 1/5*(2*(dps_A + dps_B) + dps_C)';
% ar.p(ar.qFit==1) = p_tmp_ppl + dps_C;
% % 
% [~, dps_D] = get_ppl(t_tmp, stepsize, stepsize);
% ar.p(ar.qFit==1) = p_tmp_ppl;
%dps = 1/6*(dps_A + 2*(dps_B+dps_C) + dps_D)';

% dps = dps_C';

dx = [dx; dps];

    function [fRHS_ppl, dp, sres_y] = get_ppl(t_get, step, RHS_t)
        if(~isnan(RHS_t))           
            arLink(true,t_tmp+RHS_t); 
            arChi2(false, ar.p(ar.qFit==1),1);
            [~,it] = min(abs(ar.model(m).condition(c).tExp-t_tmp-RHS_t));
            fRHS_ppl = ar.model(m).condition(c).dxdts(it,jx);
        end
        ar.p(ar.qFit==1) = p_tmp_ppl;
        [~, xSim, xSim2, xSim3, it] = PPL_chi2(t_get, true, m, c, jx, takeY, qLog10, true, step, xFit_tmp, xstd, NaN);
        if(isnan(RHS_t))
            if(takeY)
                fRHS_ppl = ar.model(m).data(c).dydt(it,jx);  
            else
                fRHS_ppl = ar.model(m).condition(c).dxdts(it,jx);  
            end
            
        end
        res = ar.res;
        sres = ar.sres;
        if(takeY)                    
           res_y = res(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it);
           sres_y = -sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,ar.qFit==1) * xstd;
           res(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it) = [];
           sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:) = [];
            grad_woy = 2*res*sres(:,ar.qFit==1);
        else
           res_y = (ar.ppl.xFit_tmp - xSim).^2 ./ xstd.^2;                     
            if(ar.config.useSensis)
                sres = ar.sres(:, ar.qFit==1);               
                if(~takeY)        
                    sxSim = zeros(1,length(ar.p));
                    sxSim(ar.model(m).condition(c).pLink) = ...
                        squeeze(ar.model(m).condition(c).sxExpSimu(it,jx,:))';
                    for j=find(ar.qLog10==1)
                        sxSim(j) = sxSim(j) * 10.^ar.p(j) * log(10);
                    end

                    if(qLog10)
                        sxSim = sxSim / 10^xSim / log(10);
                    end

                    sres_y = sxSim(ar.qFit==1);
                end
            end
            grad_woy = 2*res*sres;
        end
        
        if(takeY)
            Aorth = ar.model(m).data(c).ppl.xSens_tmp - (ar.model(m).data(c).ppl.xSens_tmp * grad_woy')./norm(grad_woy)*(grad_woy./norm(grad_woy));
        else
            Aorth = ar.model(m).condition(c).ppl.xSens_tmp - (ar.model(m).condition(c).ppl.xSens_tmp * grad_woy')./norm(grad_woy)*(grad_woy./norm(grad_woy));
        end
        Borth = sres_y - (sres_y * grad_woy')./norm(grad_woy)*(grad_woy./norm(grad_woy));
        Borth = Borth ./ norm(Borth);
        Aorth = Aorth ./ norm(Aorth);
        test_vec = Borth;
        
        if(~takeY && ((ar.model(m).condition(c).ppl.xSens_tmp * test_vec' < 0 && dir==1) || (ar.model(m).condition(c).ppl.xSens_tmp * test_vec' > 0 && dir==-1)) ...
                || takeY && ((ar.model(m).data(c).ppl.xSens_tmp * test_vec' < 0 && dir==1) || (ar.model(m).data(c).ppl.xSens_tmp * test_vec' > 0 && dir==-1)))
            test_vec = -test_vec;
        end
        
%         Borth = sres_y - (sres_y * grad_woy')./norm(grad_woy)*(grad_woy./norm(grad_woy));
%         Borth = Borth ./ norm(Borth);
        %test_vec = Borth;
        
%         Borth_perpA = (Borth - (Borth * Aorth')./norm(Aorth)*(Aorth./norm(Aorth)));
%         Aorth_perpB = -(Aorth - (Borth * Aorth')./norm(Borth)*(Borth./norm(Borth)));
%         
%         test_vec = Borth_perpA;
        %test_vec = Borth_perpA / norm(Borth_perpA) + Aorth_perpB / norm(Aorth_perpB);
        
%         if((sres_y * test_vec' > 0 && dir==1) || (sres_y * test_vec' < 0 && dir==-1))
%             test_vec = -test_vec;
%         end
        
%         if(isfield(ar.ppl,'fRHS_ppl'));
%             fRHS_ppl = ar.ppl.fRHS_ppl;
%         else
%             fRHS_ppl = ar.model(m).condition(c).dxdts(it,jx); 
%         end
%         if(jt>2)
%             curve_xFit = (ar.model(m).condition(c).ppl.x_high(jt, jx) - 2*ar.model(m).condition(c).ppl.x_high(jt-1, jx) ...
%                 + ar.model(m).condition(c).ppl.x_high(jt-2, jx)) / stepsize^2;
% 
%             if((curve_xFit>0 && dir ==1) || (curve_xFit<0 && dir==-1))
%                 if((fRHS_ppl>0 && dir==1) || (fRHS_ppl<0 && dir==-1))
%                     fRHS_ppl = (xSim - xSim3)/1.e-3;
%                     %fRHS_ppl = fRHS*1.1;
%                 else
%                     fRHS_ppl = (xSim - xSim3)/1.e-3;
%                     %fRHS_ppl = fRHS*0.9;
%                 end
%             end
%         end                

        dy_ppl = zeros(size(ar.p));
        
            if(takeY)
                %dy_ppl(ar.model(m).condition(c).pLink) = squeeze(ar.model(m).data(c).dsydt(it,jx,:))';
                dtheta_fun = matlabFunction(xSim + fRHS_ppl*stepsize - (beta.*(ar.model(m).data(c).ppl.xSens_tmp*test_vec') + xSim2));
                alpha = fzero(@(beta) dtheta_fun(beta),interval,optimoptions('fsolve','Display','off'));
            else
                %dy_ppl(ar.model(m).condition(c).pLink) = squeeze(ar.model(m).condition(c).ddxdtdps(it,jx,:))';
                dtheta_fun = matlabFunction(xSim + fRHS_ppl*stepsize - (beta.*(ar.model(m).condition(c).ppl.xSens_tmp*test_vec') + xSim2));
                % + dy_ppl(ar.qFit==1)*beta*test_vec'*stepsize
                alpha = fzero(@(beta) dtheta_fun(beta),interval,optimoptions('fsolve','Display','off'));
            end
        
        if(isnan(alpha) || ( alpha<0 && dir==1 ) || ( alpha<0 && dir==-1 ) )
            alpha=0.1/sqrt(length(ar.p(ar.qFit==1)))/gamma_tmp;
        end
        dp = alpha/abs(stepsize) .* test_vec;
    end


end