% chi2 = arPPL_intCorrection(general_struct, t_tmp, dir)
% 
%   This function is used by arIntegratePredBand.m
% 
% general_struct    PPL result struct
% t_tmp             The time point of the prediction, passed later to arLink
% dir               direction, either 1 or -1
% 
% 
% Written by Helge, tried to be documented by Clemens.

function chi2 = arPPL_intCorrection(general_struct, t_tmp, dir)
    global ar;
    
    %fill temporary variables
    data_cond = general_struct.data_cond;
    x_y = general_struct.x_y;
    m=general_struct.m;
    c=general_struct.c; 
    jx=general_struct.jx;     
    xstd = ar.ppl.options.xstd;
           
    chi2_best = ar.ppl.chi2_threshold - ar.ppl.dchi2;
    bkp_fit = ar.ppl.xFit_tmp;
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
    x_orig = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);
    
    chi2_orig=arPPL_Chi2Corr(general_struct, t_tmp);
    
    fit_a = (chi2_orig-chi2_best) / (bkp_fit-x_orig)^2;
    if(fit_a>0)
        ar.ppl.xFit_tmp = x_orig + dir*sqrt(ar.ppl.dchi2/fit_a);
    else
        ar.ppl.xFit_tmp = x_orig + dir*xstd;
    end
    fit_quadr = ar.ppl.xFit_tmp;
    chi2=arPPL_Chi2Corr(general_struct, t_tmp);
    chi2_quadr = chi2;
    
    %Fit second order polynomial
    if((chi2 - ar.ppl.chi2_threshold) > 0.2 || (chi2 - ar.ppl.chi2_threshold) < -0.2 )
        x_fit = [0; ar.ppl.xFit_tmp-x_orig; bkp_fit-x_orig];
        y_fit = [chi2_best; chi2_quadr; chi2_orig];
        cf = fit(x_fit, y_fit, 'poly2', 'Lower', [0 -Inf chi2_best*0.95], 'Upper', [Inf Inf chi2_best*1.05], 'Robust', 'off');
        if(abs(cf.p1)>1.e-3)
            ar.ppl.xFit_tmp=x_orig -cf.p2/(2*cf.p1) + dir*(sqrt((cf.p2/(2*cf.p1))^2 - (cf.p3-ar.ppl.chi2_threshold)/cf.p1));
        else
            ar.ppl.xFit_tmp = x_orig + (ar.ppl.chi2_threshold - cf.p3)/cf.p2;
        end
    	chi2=arPPL_Chi2Corr(general_struct, t_tmp);
        chi2_poly = chi2;
        fit_poly = ar.ppl.xFit_tmp;
    end
    
    
    %Fit cubic spline
    if((chi2 - ar.ppl.chi2_threshold) > 0.2 || (chi2 - ar.ppl.chi2_threshold) < -0.2 )    
        x_fit = [0; bkp_fit-x_orig; fit_poly-x_orig; fit_quadr-x_orig];
        y_fit = [chi2_best; chi2_orig; chi2_poly; chi2_quadr];
        ar.ppl.xFit_tmp = interp1(x_fit,y_fit,chi2_best+ar.ppl.dchi2,'pchip','extrap');
        chi2=arPPL_Chi2Corr(general_struct, t_tmp);
    end
end