function chi2 = PPL_doSim_calc(chi2_orig, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL)
    global ar;
    chi2_best = ar.ppl.chi2_95 - ar.ppl.dchi2 - 0.5;
    bkp_fit = ar.ppl.xFit_tmp;
    bkp_p = ar.p(ar.qFit==1);
    %chi2_orig
    if(doPPL)
        x_corr = 4;
    else
        x_corr = 0;
    end
    chi2_orig=PPL_corr(x_corr, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2_orig);
    
    fit_a = (chi2_orig-chi2_best) / (bkp_fit-x_orig)^2;
    if(fit_a>0)
        ar.ppl.xFit_tmp = x_orig + dir*sqrt(ar.ppl.dchi2/fit_a);
    else
        ar.ppl.xFit_tmp = x_orig + dir*xstd;
    end

    chi2=PPL_corr(x_corr, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2_orig);
    
    if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
        x_fit = [0; ar.ppl.xFit_tmp-x_orig; bkp_fit-x_orig];
        y_fit = [chi2_best; chi2; chi2_orig];
        cf = fit(x_fit, y_fit, 'poly2', 'Lower', [0 -Inf chi2_best*0.95], 'Upper', [Inf Inf chi2_best*1.05], 'Robust', 'off');
        if(cf.p1~=0)
            ar.ppl.xFit_tmp=x_orig -cf.p2/(2*cf.p1) + dir*(sqrt((cf.p2/(2*cf.p1))^2 - (cf.p3-ar.ppl.chi2_95 + 0.5)/cf.p1));
        else
            ar.ppl.xFit_tmp = x_orig + (ar.ppl.chi2_95 - 0.5 - cf.p3)/cf.p2;
        end
        %ar.ppl.xFit_tmp
        chi2=PPL_corr(x_corr, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);
    
        %chi2
        if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
            x_fit = [x_fit; ar.ppl.xFit_tmp-x_orig];
            y_fit = [y_fit; chi2];
            ft_ = fittype('A*X^3+B*X^2+C*X+D',...
                 'dependent',{'y'},'independent',{'X'},...
                 'coefficients',{'A', 'B', 'C', 'D'});
            cf = fit(x_fit, y_fit, ft_, 'Lower', [-Inf -Inf -Inf chi2_best*0.95], 'Upper', [Inf Inf Inf chi2_best*1.05],'StartPoint',[0 1 0 chi2_best], 'Robust', 'off');
            cf.D = cf.D -ar.ppl.chi2_95 + 0.5 ;
            find_y = fzero(cf,x_fit(end));
            ar.ppl.xFit_tmp = x_orig + find_y;
            chi2=PPL_corr(x_corr, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);
    
            if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                x_fit = [x_fit; ar.ppl.xFit_tmp-x_orig];
                y_fit = [y_fit; chi2];
                [sort_y, ind_y] = sort(y_fit);
                ar.ppl.xFit_tmp = x_orig + dir*feval(cf,ar.ppl.chi2_95-0.5);
                
                chi2=PPL_corr(x_corr, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);
            end
        end
    end
end