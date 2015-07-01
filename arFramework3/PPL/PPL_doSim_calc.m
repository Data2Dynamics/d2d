function chi2 = PPL_doSim_calc(chi2_orig, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL)
    global ar;
    chi2_best = ar.ppl.chi2_95 - ar.ppl.dchi2 - 0.5;
    bkp_fit = ar.ppl.xFit_tmp;
    bkp_p = ar.p(ar.qFit==1);
    %chi2_orig
    if(doPPL)
        [chi2_orig,~,exitflag]=PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2_orig);
    else
        chi2_orig=PPL_corr(0, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, 0);
    end
    %ar.ppl.xFit_tmp
    %chi2_orig
%     x_fit = [0; bkp_fit-x_orig];
%     y_fit = [0; chi2_orig-chi2_best];
%     ft_ = fittype({'x^2','x'},...
%         'coefficients',{'p1','p2'});
%     cf = fit(x_fit, y_fit, ft_, 'Lower', [0 -Inf], 'Upper', [Inf Inf], 'Robust', 'off');        
%     if(cf.p1~=0)
%         ar.ppl.xFit_tmp=x_orig -cf.p2/(2*cf.p1) + dir*(sqrt((cf.p2/(2*cf.p1))^2  + ar.ppl.dchi2/cf.p1));
%     else
%         ar.ppl.xFit_tmp = x_orig + ar.ppl.dchi2/cf.p2;
%     end
    fit_a = (chi2_orig-chi2_best) / (bkp_fit-x_orig)^2;
    if(fit_a>0)
        ar.ppl.xFit_tmp = x_orig + dir*sqrt(ar.ppl.dchi2/fit_a);
    else
        ar.ppl.xFit_tmp = x_orig + dir*xstd;
    end
    %ar.ppl.xFit_tmp
    %fit_a = (chi2_orig-chi2_best) / (bkp_fit-x_orig);
    %ar.ppl.xFit_tmp = x_orig + ar.ppl.dchi2/fit_a;
    if(doPPL)
        [chi2,~,exitflag]=PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2_orig);
    else
        chi2=PPL_corr(0, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, 0);
    end
    %chi2
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
        if(doPPL)
            chi2=PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);
        else
            chi2=PPL_corr(0, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, 0);
        end
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
            if(doPPL)
                chi2=PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);
            else
                chi2=PPL_corr(0, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, 0);
            end
            %ar.ppl.xFit_tmp
            %chi2
            if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                x_fit = [x_fit; ar.ppl.xFit_tmp-x_orig];
                y_fit = [y_fit; chi2];
                [sort_y, ind_y] = sort(y_fit);
                %cf = fit(x_fit, y_fit, ft_,fo_);
%                cf = fit(sort_y, x_fit(ind_y), 'cubicinterp');%, 'Lower', [0 -Inf chi2_best*0.95], 'Upper', [Inf Inf chi2_best*1.05], 'Robust', 'off');
    %             if(cf.p1~=0)
    %                 ar.ppl.xFit_tmp=x_orig -cf.p2/(2*cf.p1) + dir*(sqrt((cf.p2/(2*cf.p1))^2 - (cf.p3-ar.ppl.chi2_95 + 0.5)/cf.p1));
    %             else
    %                 ar.ppl.xFit_tmp = x_orig + (ar.ppl.chi2_95 - 0.5 - cf.p3)/cf.p2;
    %             end
                ar.ppl.xFit_tmp = x_orig + dir*feval(cf,ar.ppl.chi2_95-0.5);
                
                %cf = fit(x_fit, y_fit, ft_,fo_);
                
%                 if(cf.p1~=0)
%                     ar.ppl.xFit_tmp=x_orig -cf.p2/(2*cf.p1) + dir*(sqrt((cf.p2/(2*cf.p1))^2 - (cf.p3-ar.ppl.chi2_95 + 0.5)/cf.p1));
%                 else
%                     ar.ppl.xFit_tmp = x_orig + (ar.ppl.chi2_95 - 0.5 - cf.p3)/cf.p2;
%                 end
                
                if(doPPL)
                    chi2=PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);
                else
                    chi2=PPL_corr(0, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, 0);
                end
                %chi2
            end
        end
    end
end