function chi2 = PPL_doSim_calc(m, c, jx, xstd, t_tmp, qLog10, dir, takeY)
    global ar;
    if(takeY)
        data_cond = 'data';
        x_y = 'y';
    else
        data_cond = 'condition';
        x_y = 'x';
    end
    
    chi2_best = ar.ppl.chi2_95 - ar.ppl.dchi2 - 0.5;
    bkp_fit = ar.ppl.xFit_tmp;
    bkp_p = ar.p(ar.qFit==1);
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
    x_orig = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);
    
    chi2_orig=PPL_corr(m, c, jx, xstd, t_tmp, qLog10, takeY);
    
    fit_a = (chi2_orig-chi2_best) / (bkp_fit-x_orig)^2;
    if(fit_a>0)
        ar.ppl.xFit_tmp = x_orig + dir*sqrt(ar.ppl.dchi2/fit_a);
    else
        ar.ppl.xFit_tmp = x_orig + dir*xstd;
    end
    fit_quadr = ar.ppl.xFit_tmp;
    chi2=PPL_corr(m, c, jx, xstd, t_tmp, qLog10, takeY);
    chi2_quadr = chi2;
    
    %Fit second order polynomial
    if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
        x_fit = [0; ar.ppl.xFit_tmp-x_orig; bkp_fit-x_orig];
        y_fit = [chi2_best; chi2_quadr; chi2_orig];
        cf = fit(x_fit, y_fit, 'poly2', 'Lower', [0 -Inf chi2_best*0.95], 'Upper', [Inf Inf chi2_best*1.05], 'Robust', 'off');
        if(cf.p1~=0)
            ar.ppl.xFit_tmp=x_orig -cf.p2/(2*cf.p1) + dir*(sqrt((cf.p2/(2*cf.p1))^2 - (cf.p3-ar.ppl.chi2_95 + 0.5)/cf.p1));
        else
            ar.ppl.xFit_tmp = x_orig + (ar.ppl.chi2_95 - 0.5 - cf.p3)/cf.p2;
        end
    	chi2=PPL_corr(m, c, jx, xstd, t_tmp, qLog10, takeY);
        chi2_poly = chi2;
        fit_poly = ar.ppl.xFit_tmp;
    end
    
    
    %Fit cubic spline
    if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )    
        x_fit = [0; bkp_fit-x_orig; fit_poly-x_orig; fit_quadr-x_orig];
        y_fit = [chi2_best; chi2_orig; chi2_poly; chi2_quadr];
        ar.ppl.xFit_tmp = interp1(x_fit,y_fit,chi2_best+ar.ppl.dchi2,'pchip','extrap');
        chi2=PPL_corr(m, c, jx, xstd, t_tmp, qLog10, takeY);
    end
end