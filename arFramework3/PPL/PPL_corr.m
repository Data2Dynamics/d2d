function [chi2_out, xSim, exitflag] = PPL_corr(x, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, x_factor)
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
            if(dir==1)
                lb_corr=[ar.lb(ar.qFit==1) x_orig];
                ub_corr=[ar.ub(ar.qFit==1) Inf];
            else
                lb_corr=[ar.lb(ar.qFit==1) -Inf];
                ub_corr=[ar.ub(ar.qFit==1) x_orig];
            end
            if(x == 3)
                constr_option=optimoptions('fmincon','Algorithm','interior-point','Display','none','GradConstr','on','GradObj','on','TolCon',1.e-5,'MaxIter',50,'TolFun',0,'TolX',1.e-8,'Hessian','user-supplied','HessFcn',@Corr_Hess,'InitTrustRegionRadius',sqrt(0.1*length(ar.p(ar.qFit==1))),'AlwaysHonorConstraints','none');                
                [pFit, chi2, exitflag] = ...
                    fmincon(@ppl_constr_fkt, [ar.p(ar.qFit==1) ar.ppl.xFit_tmp], [], [],[], [], lb_corr, ub_corr, @constr_lsq, constr_option, x);   
                %exitflag
            elseif(x == 4)
                
                if(doPPL)
                    constr_option2=optimoptions('fmincon','Algorithm','interior-point','Display','none','GradConstr','on','GradObj','on','TolCon',1.e-5,'MaxIter',60,'TolFun',0,'TolX',1.e-8,'Hessian','user-supplied','HessFcn',@Corr_Hess,'InitTrustRegionRadius',sqrt(0.1*length(ar.p(ar.qFit==1))),'AlwaysHonorConstraints','none','SubproblemAlgorithm','cg'); %
                    [pFit, chi2, exitflag] = ...
                        fmincon(@ppl_constr_fkt, ar.p(ar.qFit==1), [], [],[], [], ar.lb(ar.qFit==1), ar.ub(ar.qFit==1), @constr_lsq, constr_option2, x);
                else
                    corr_option=optimoptions('fmincon','Algorithm','trust-region-reflective','DerivativeCheck','off','Display','none','GradObj','on','TolCon',1.e-2,'MaxIter',50,'TolFun',1.e-5,'TolX',1.e-4,'Hessian','user-supplied');%,'OutputFcn',@ppl_corr);
                    [pFit, chi2, exitflag] = ...
                        fmincon(@ppl_corr_fkt, ar.p(ar.qFit==1), [], [],[], [], ar.lb(ar.qFit==1), ar.ub(ar.qFit==1), [], corr_option, x);                  
                end
                
                  
            elseif(x==0)
                    corr_option=optimoptions('fmincon','Algorithm','trust-region-reflective','DerivativeCheck','off','Display','none','GradObj','on','TolCon',1.e-2,'MaxIter',50,'TolFun',0,'TolX',1.e-4,'Hessian','user-supplied');
                    [pFit, chi2, exitflag] = ...
                        fmincon(@ppl_corr_fkt, ar.p(ar.qFit==1), [], [],[], [], ar.lb(ar.qFit==1), ar.ub(ar.qFit==1), [], corr_option, x);   
            elseif(x==1)
                corr_option=optimoptions('fmincon','Algorithm','trust-region-reflective','DerivativeCheck','off','Display','none','GradObj','on','TolCon',1.e-2,'MaxIter',50,'TolFun',0,'TolX',1.e-4,'Hessian','user-supplied');
                [pFit, chi2, exitflag] = ...
                    fmincon(@ppl_corr_fkt, ar.ppl.xFit_tmp, [], [],[], [], lb_corr(end), ub_corr(end), [], corr_option, x);   
            elseif(x==2)
                corr_option=optimoptions('fmincon','Algorithm','trust-region-reflective','DerivativeCheck','off','Display','none','GradObj','on','TolCon',1.e-2,'MaxIter',50,'TolFun',1.e-6,'TolX',1.e-6,'Hessian','user-supplied');%,'OutputFcn',@ppl_corr);
                [pFit, chi2, exitflag] = ...
                    fmincon(@ppl_corr_fkt, [ar.p(ar.qFit==1) ar.ppl.xFit_tmp], [], [],[], [], lb_corr, ub_corr, [], corr_option, x);                  
            end
            if(x==2 || x == 3)
                ar.p(ar.qFit==1) = pFit(1:end-1);
                if(x==3) 
                    ar.ppl.xFit_tmp = xSim;
                    chi2_out = chi2-0.75;
                    %ar.ppl.grad_tmp = chi2_grad;
                else
                    ar.ppl.xFit_tmp = pFit(end);
                    chi2_out = chi2-0.25;
                end
                
            elseif(x==1)
                ar.ppl.xFit_tmp = pFit;
                chi2_out = chi2-0.25;
                %ar.ppl.xFit_tmp
                
            elseif(x==0 || x == 4 || x==5)
                ar.p(ar.qFit==1) = pFit;   
                chi2_out = chi2;
                %ar.ppl.grad_tmp = chi2_grad;
            end
            if(exitflag < 1)
                outputstr = '';
                switch exitflag                   
                    case 0
                        outputstr = 'Too many function evaluations or iterations.';
                    case -1
                        outputstr = 'Stopped by output/plot function.';
                    case -2
                        outputstr = 'Bounds are inconsistent.';
                    case -3
                        outputstr = 'Objective function below limit with constraints not satisfied.';                   
                end
                %fprintf('fmincon finished after %i interations: %s\n', ar.fit.output.iterations, outputstr);
            end

    function [chi2, grad, H] = ppl_corr_fkt(pTrial, x)

        iter = iter+1;
        if(x==2)
            xExp_tmp=pTrial(end);
            pTrial=pTrial(1:end-1);
        elseif(x==1)
            xExp_tmp=pTrial;
            pTrial = ar.p(ar.qFit==1);
        elseif(x==0 || x==4)
            xExp_tmp=ar.ppl.xFit_tmp;
        end
        if (x==1 || ~isequal(pTrial,Last_pars))
            arLink(true,t_tmp,takeY,jx, c, m,xExp_tmp,xstd);  
            try
                arChi2(ar.config.useSensis, pTrial,1)
            catch
                fprintf('arChi2 failed at t= %d and xFit = %d for x=%i, reset xFit! \n', t_tmp, ar.ppl.xFit_tmp, x);
                if(dir==1)
                    ar.ppl.xFit_tmp = ar.model(m).(data_cond)(c).ppl.x_high_vpl(jt, jx);                                       
                else
                    ar.ppl.xFit_tmp = ar.model(m).(data_cond)(c).ppl.x_low_vpl(jt, jx);    
                end
                xExp_tmp = ar.ppl.xFit_tmp;
                arLink(true,t_tmp,true,jx, c, m,xExp_tmp,xstd);  
                arChi2(ar.config.useSensis, pTrial,1);
            end

            %res = [ar.res ar.constr];
            res = ar.res;
            
           [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
           if(takeY && length(find(ar.model(m).data(c).tExp==t_tmp))>1)
                it = it+1;               
           end
           xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx); 
            
            if(qLog10)
                xSim = log10(xSim);
            end
            if(nargout>1 && ar.config.useSensis)
                sres = ar.sres(:, ar.qFit==1);
            end

            if(~takeY)
    %              if(~doPPL)          
                res(end+1) = (xExp_tmp-xSim)/xstd;   
                if(x~=1)
                    ar.res(end+1) = (xExp_tmp-xSim)/xstd; 
                end
                if(nargout>1 && ar.config.useSensis)% && ~doPPL)

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
                    if(x~=1)
                        ar.sres(end+1,:) = -sxSim / xstd;
                    end
                end
            end
            Last_pars = pTrial;
        else
            res = ar.res;
            sres = ar.sres(:,ar.qFit==1);          
        end
        if(nargout>1 && ar.config.useSensis)
            if(takeY)
               sres_tmp = ar.sres(ar.ppl.resi_tmp,ar.qFit==1);          
            else
                sres_tmp = sres(end,:);
            end
           
            if(x~=1)
                grad = 2*res*sres;
                grad_tmp = grad;
                H = 2*(sres') * sres;
            elseif(x==1)
                grad = 2*(xExp_tmp-xSim)/xstd^2;
                grad_tmp = grad;
                H = 2./(xstd.^2);
            end
            if(x==2 && mod(iter,2)==0)
                grad_tmp = [grad, 2*(xExp_tmp-xSim)/xstd^2];
                H_hinten=2*nansum(sres_tmp,1)'./xstd;
                H_unten=[2*nansum(sres_tmp,1)./xstd, 2./(xstd.^2)];

                H= [H, H_hinten];
                H= [H; H_unten];
            elseif(x==2)
                grad = [grad, 0];
                H = [H, zeros(size(H,1),1)];
                H = [H; zeros(1,size(H,2))];
            end
            if(x==1 || (x==2 && mod(iter,2)==0))                
                grad = (1+2*(nansum(res.^2)-ar.ppl.chi2_95))*grad_tmp;   
                H = (1+2*(nansum(res.^2)-ar.ppl.chi2_95))*H + 2*(grad_tmp')*grad_tmp;
            end
            %grad = (1+2*(ar.chi2fit-ar.ppl.chi2_95))*grad_tmp;            
            %H = (1+2*(ar.chi2fit-ar.ppl.chi2_95))*H + 2*(grad_tmp')*grad_tmp;
        end        
        chi2 = nansum(res.^2);
        if(x==1 || (x==2 && mod(iter,2)==0))  
            chi2 = nansum(res.^2) + (nansum(res.^2)-ar.ppl.chi2_95)^2; 
        end
        %chi2 = ar.chi2fit + (ar.chi2fit-ar.ppl.chi2_95)^2;    
    end

    function [chi2, grad] = ppl_constr_fkt(pTrial, x)
    iter = iter+1;
    if(x == 3)
        xExp_tmp=pTrial(end);
        pTrial=pTrial(1:end-1);
    else
        xExp_tmp = ar.ppl.xFit_tmp;
    end
    if ~isequal(pTrial,Last_pars)
        arLink(true,t_tmp,takeY,jx, c, m,xExp_tmp,xstd);  
        try
            arChi2(ar.config.useSensis, pTrial,1);
        catch
            fprintf('arChi2 failed at t= %d and xFit = %d, reset xFit! \n', t_tmp, ar.ppl.xFit_tmp);
            if(~doPPL)
                if(dir==1)
                    ar.ppl.xFit_tmp = ar.model(m).(data_cond)(c).ppl.x_high_vpl(jt, jx);   
                else
                    ar.ppl.xFit_tmp = ar.model(m).(data_cond)(c).ppl.x_low_vpl(jt, jx);                                       
                end
                xExp_tmp = ar.ppl.xFit_tmp;
                arLink(true,t_tmp,takeY,jx, c, m,xExp_tmp,xstd);  
                arChi2(ar.config.useSensis, pTrial,1);
            end
        end
        %res = [ar.res ar.constr];
        res = ar.res;
        if(nargout>1 && ar.config.useSensis)
            sres = ar.sres(:, ar.qFit==1);
        end      
        [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
       if(takeY && length(find(ar.model(m).(data_cond)(c).tExp==t_tmp))>1)
            it = it+1;               
       end
        xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);
    
        if(qLog10)
            xSim = log10(xSim);
        end
        if(nargout>1 && ar.config.useSensis)

            if(~takeY)
                res(end+1) = (xExp_tmp-xSim)/xstd; 

                ar.res(end+1) = (xExp_tmp-xSim)/xstd;  
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
                    ar.sres(end+1,:) = -sxSim / xstd;
                end                
            end            

        end
        
        if(x==4 && takeY)
           sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:) = []; 
           res(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it) = [];   
        elseif(x==4 && ~takeY)
            res(end) = [];
            sres(end,:) = [];
        end
        
        grad = 2*res*sres;
        if(x == 3)% && mod(iter,2))
            grad = [grad, 2*(xExp_tmp-xSim)/xstd^2]; 
        end    
        
        Last_pars = pTrial;
     else
        res = ar.res;
        sres = ar.sres(:,ar.qFit==1);
        if(x==4 && takeY)
           sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:) = []; 
           res(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it) = [];   
        elseif(x==4 && ~takeY)
            res(end) = [];
            sres(end,:) = [];
        end
        if(x == 3)% && mod(iter,2))
            grad = [2*res*sres, 2*(xExp_tmp-xSim)/xstd^2]; 
        else
            grad = 2*res*sres;
        end
     end
        chi2 = nansum(res.^2);    
        if(doPPL && x ==3)%  && mod(iter,2))
            chi2 = chi2 + (nansum(res.^2)-ar.ppl.chi2_95-0.5)^2; 
            grad = (1+2*(nansum(res.^2)-ar.ppl.chi2_95-0.5))*grad;  
        elseif(x  == 3)% && mod(iter,2))
            chi2 = chi2 + (nansum(res.^2)-ar.ppl.chi2_95)^2;
            grad = (1+2*(nansum(res.^2)-ar.ppl.chi2_95))*grad;  
        end                                
    end

    function [H] = Corr_Hess(xs,lambda, x)
        if ((~isequal(xs(1:end-1),Last_pars) && x==3) || (~isequal(xs,Last_pars) && x==4))
            if(takeY && x==3)
               arLink(true,t_tmp,true,jx, c, m,xs(end),xstd);  
            elseif(takeY && x==4)
                arLink(true,t_tmp,true,jx, c, m,ar.ppl.xFit_tmp,xstd);
            else
                arLink(true,t_tmp);
            end
            if(x == 3)
                arChi2(ar.config.useSensis, xs(1:end-1),1);
            else
                arChi2(ar.config.useSensis, xs,1);
            end
            res = ar.res;
            [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
           if(takeY && length(find(ar.model(m).(data_cond)(c).tExp==t_tmp))>1)
                it = it+1;               
           end
            xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);
           
            if(qLog10)
                xSim = log10(xSim);
            end
            sres = ar.sres(:, ar.qFit==1);         
    %             if(~isempty(ar.sconstr))
    %                 sres = [ar.sres(:, ar.qFit==1); ar.sconstr(:, ar.qFit==1)];
    %             end
            if(~takeY)
                res(end+1) = (xs(end)-xSim)/xstd;  
                ar.res(end+1) = (xs(end)-xSim)/xstd;  
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
                ar.sres(end+1,:) = -sxSim / xstd;
            end
            if(x==3)
                Last_pars = xs(1:end-1);
            else
                Last_pars = xs;
            end
        else
            sres = ar.sres(:, ar.qFit==1);
            res = ar.res;
        end
        if(takeY)
           sres_tmp = sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:); 
        else
            sres_tmp = sres(end,:);
        end
        ar.ppl.sres_tmp = sres_tmp;
        if(x == 4)
            if(takeY)
                res(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it) = [];
                sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:) = [];
            else
                res(end) = [];
                sres(end,:) = [];               
            end
        end
        grad = 2*res*sres;
        H = 2*(sres)' * sres;
        if(x == 3)% && mod(iter,2))        
            grad_tmp = [grad, 2*(xs(end)-xSim)/xstd^2];

            H_hinten=2*nansum(sres_tmp,1)'./xstd;
            H_unten=[2*nansum(sres_tmp,1)./xstd, 2./(xstd.^2)];

            H= [H, H_hinten];
            H= [H; H_unten];

            if(doPPL)
                H = (1+2*(nansum(res.^2)-ar.ppl.chi2_95-0.5))*H + 2*(grad_tmp')*grad_tmp;
            else
                H = (1+2*(nansum(res.^2)-ar.ppl.chi2_95))*H + 2*(grad_tmp')*grad_tmp;
            end

            H_eq = 2*sres_tmp'*sres_tmp;
            H_eq_hinten=2*nansum(sres_tmp,1)'./xstd;
            H_eq_unten=[2*nansum(sres_tmp,1)./xstd, 2./(xstd.^2)];
            H_eq = [H_eq, H_eq_hinten];
            H_eq = [H_eq; H_eq_unten];  
            H=H+lambda.eqnonlin*H_eq;
        else
            H=H+lambda.eqnonlin*2*(sres_tmp'*sres_tmp);
        end
        
    end
    function [cineq,ceq,GC,GCeq] = constr_lsq(p, x)
        if ((~isequal(p(1:end-1),Last_pars) && x==3) || (~isequal(p,Last_pars) && x==4))
            if(takeY && x==3)
               arLink(true,t_tmp,true,jx, c, m,p(end),xstd);  
            elseif(takeY && x==4)
                arLink(true,t_tmp,true,jx, c, m,ar.ppl.xFit_tmp,xstd);
            else
                arLink(true,t_tmp);
            end
            if(x == 3)
                arChi2(ar.config.useSensis, p(1:end-1),1); 
            elseif(x == 4)
                arChi2(ar.config.useSensis, p,1); 
            end
            %res = [ar.res ar.constr];   
            res = ar.res;
            [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
           if(takeY && length(find(ar.model(m).(data_cond)(c).tExp==t_tmp))>1)
                it = it+1;               
           end
            xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);          
           
            if(qLog10)
                xSim = log10(xSim);
            end
            sres = ar.sres(:, ar.qFit==1);         
            if(~takeY)
                res(end+1) = (p(end)-xSim)/xstd;  
                ar.res(end+1) = (p(end)-xSim)/xstd;  
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
                ar.sres(end+1,:) = -sxSim / xstd;
            end
            if(x==3)
                Last_pars = p(1:end-1);
            else
                Last_pars = p;
            end
        else
            res = ar.res;
            sres = ar.sres(:, ar.qFit==1);
        end
        
        cineq=[];
        GC=[];    
        if(takeY)
           sres_tmp = sres(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it,:);   
           res_tmp = res(ar.res_mLink==m & ar.res_type==1 & ar.res_dLink==c & ar.res_yLink==jx & ar.res_tLink==it);
        else
            sres_tmp = sres(end,:);
            res_tmp = res(end);
        end
        %ceq=ar.chi2fit+(p(end)-xSim)^2/xstd^2-ar.ppl.chi2_95;
        if(x==3)
            ceq= (xSim-p(end))^2./xstd^2 - 0.5;
            grad = 2*res_tmp*sres_tmp;
            grad = [grad, 2*(p(end)-xSim)./xstd.^2];

            GCeq=grad';
        else
            ceq= (xSim-ar.ppl.xFit_tmp)^2/xstd^2;
            GCeq = (-2*(xSim-ar.ppl.xFit_tmp)*sres_tmp/xstd)';
        end
    end
end