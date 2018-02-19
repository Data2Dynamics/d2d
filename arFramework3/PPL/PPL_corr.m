function [chi2_out, xSim, exitflag] = PPL_corr(m, c, jx, xstd, t_tmp, qLog10, takeY)
        global ar;
        xSim = NaN;
        it = NaN;
        Last_pars = [];
        iter = 0;
        if(takeY)
            data_cond = 'data';
            x_y = 'y';
        else
            data_cond = 'condition';
            x_y = 'x';
        end
        %try
        xExp_tmp=ar.ppl.xFit_tmp;
        arLink(true,t_tmp,takeY,jx, c, m,xExp_tmp,xstd);  
        
    
        [pFit, ~, ~, exitflag] = ...
            lsqnonlin(@ppl_merit_fkt, ar.p(ar.qFit==1), ar.lb(ar.qFit==1), ar.ub(ar.qFit==1), ar.config.optim);
        ar.p(ar.qFit==1) = pFit;
        if(exitflag < 1)
            outputstr = '';
            switch exitflag
                case 1
                    outputstr = 'LSQNONLIN converged to a solution.';
                case 2
                    outputstr = 'Change in X too small.';
                case 3
                    outputstr = 'Change in RESNORM too small.';
                case 4
                    outputstr = 'Computed search direction too small.';
                case 0
                    outputstr = 'Too many function evaluations or iterations.';
                case -1
                    outputstr = 'Stopped by output/plot function.';
                case -2
                    outputstr = 'Bounds are inconsistent.';
                case -3
                    outputstr = 'Regularization parameter too large (Levenberg-Marquardt).';
                case -4
                    outputstr = 'Line search failed.';
            end
            %fprintf('lsqnonlin finished after %i interations: %s\n', ar.fit.output.iterations, outputstr);
        end 
        arCalcMerit(0, ar.p(ar.qFit==1),1)     
        if(~takeY)
            chi2_out = arGetMerit('chi2') + ((xExp_tmp-xSim)/xstd).^2;
        else
            chi2_out = arGetMerit('chi2');
        end

        function [res, sres] = ppl_merit_fkt(pTrial)
            arCalcMerit(ar.config.useSensis, pTrial, 1)
            res = [ar.res ar.constr];
            %res = ar.res;
            if(nargout>1 && ar.config.useSensis)
                sres = ar.sres(:, ar.qFit==1);
                if(~isempty(ar.sconstr))
                    sres = [ar.sres(:, ar.qFit==1); ar.sconstr(:, ar.qFit==1)];
                end
            end
            [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
            if(takeY && length(find(ar.model(m).data(c).tExp==t_tmp))>1)
                it = it+1;               
            end
            xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx); 
            
            if(qLog10)
                xSim = log10(xSim);
            end
            if(~takeY)
                if(qLog10)
                        xSim = log10(xSim);
                end	
                res(end+1) = (xExp_tmp-xSim)/xstd;

                if(nargout>1 && ar.config.useSensis)

                    sxSim = zeros(1,length(ar.p));           
                    sxSim(ar.model(m).condition(c).pLink) = ...
                            squeeze(ar.model(m).condition(c).sxExpSimu(it,jx,:))';
                    for j10=find(ar.qLog10==1)
                        sxSim(j10) = sxSim(j10) * 10.^ar.p(j10) * log(10);
                    end
                    if(qLog10)
                        sxSim = sxSim / 10^xSim / log(10);
                    end

                    sres(end+1,:) = -sxSim(ar.qFit==1) / xstd;
                end    
            end            
        end   
   
end