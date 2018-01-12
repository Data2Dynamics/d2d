function ppl_calc(m, c, jx, xFit, p, t, doPPL, takeY, dir, stepsize, xstd, ed_steps, pReset, chi2start, backward, fineInt)
    global ar
    ar.ppl.xFit_tmp = xFit;
    chi2 = chi2start + ar.ppl.dchi2;
    xSim = NaN;
    npre=0;
    if(takeY)
        data_cond = 'data';
        x_y = 'y';
    else
        data_cond = 'condition';
        x_y = 'x';
    end
          
    if(dir==1)
        high_low = '_high';
    else
        high_low = '_low';
    end
    
    if(isnan(xFit))
        error('Starting value is NaN, check calculation of profiles! \n')
    end
    t_dir = 1;
    t_tmp=t;
    jt = 0;
    if(backward)
        jt = ar.ppl.nsteps-1;
        t_dir = -1;
    end
    qLog10 = ar.ppl.qLog10;
    if(takeY)
        xLabel = myNameTrafo(ar.model(m).data(c).y{jx});
        gamma_tmp = ar.model(m).data(c).ppl.gamma(jx);
    else
        xLabel = myNameTrafo(ar.model(m).x{jx});
        gamma_tmp = ar.model(m).condition(c).ppl.gamma(jx);
    end    
    nsteps=ar.ppl.nsteps;
    ar.p=squeeze(p)';
    if(size(ar.p,1)~=1)
        ar.p=ar.p';
    end
    tcount = 1;
    arWaitbar(0);
    while jt<nsteps
        jt = jt+1;
        corr_tmp = 0;

        [~,it_orig] = min(abs(ar.model(m).(data_cond)(c).tFine-t_tmp-t_dir*stepsize));
        
        if(toc>tcount)        
            if(dir==1)
                string_tmp = 'upper';
            else
                string_tmp = 'lower';            
            end
            arWaitbar((jt+npre), nsteps, sprintf(['PPL-integration (' string_tmp ' bound) for %s at t=%g %i/%i'], xLabel, t_tmp, jt, nsteps));
            tcount = tcount + 0.5; % update every half second
        end        

        if(ed_steps==true)              
            %VPL part    
            [chi2, ~] = PPL_chi2(t_tmp, false, m, c, jx, takeY, qLog10, t_dir*stepsize, ar.ppl.xFit_tmp, xstd);        
            [dps, dx, gamma_tmp] = VPL_Int(t_tmp, m, c, jx, takeY, qLog10, t_dir*stepsize, ar.ppl.xFit_tmp, xstd, gamma_tmp, fineInt, t_dir);             
            
            ar.p(ar.qFit==1)=ar.p(ar.qFit==1) + dps'*stepsize;
            ar.ppl.xFit_tmp = ar.ppl.xFit_tmp + dx(1)*stepsize;           
            t_tmp = t_tmp + t_dir*stepsize;              
            
            if(fineInt || mod(jt/(floor(nsteps/10)),1)==0)
                [chi2, xSim] = PPL_corr(m, c, jx, xstd, t_tmp, qLog10, takeY);                
            else
                [chi2, xSim] = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, stepsize, ar.ppl.xFit_tmp, xstd);
            end
                       
            if(sum(isnan(dx))>0)               
               ar.p=pReset;
               if(takeY)
                    arLink(true,ar.model(m).data(c).tExp(1),true,jx, c, m,NaN);
               end
               fprintf('ERROR IN STEP AT T=%d \n', t_tmp);
               return;
            end            
        else

            xSim = getxFit(ar.ppl.xFit_tmp);
            t_tmp = t_tmp + t_dir*stepsize;
            ar.ppl.xFit_tmp = xSim;
            chi2 = PPL_doSim_calc(m, c, jx, xstd, t_tmp, qLog10, dir, takeY); 

        end

        if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2)

            i_count=0;
            while((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )  
                fprintf('Correction at t=%d \n',t_tmp);
                chi2 = PPL_doSim_calc(m, c, jx, xstd, t_tmp, qLog10, dir, takeY);

                i_count = i_count + 1 ;
                if(i_count>1)
                    break;
                end
            end

            if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                fprintf('Try to restart at t=%d diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_95 + 0.5));
                ar.p = pReset;
                arLink(true,0.,takeY,jx, c, m,NaN);
                [ar.ppl.xFit_tmp, ar.p] = xstart_ppl(m, c, jx, t_tmp, doPPL, xstd, pReset, chi2start, 10, takeY, false, dir, ar.ppl.xFit_tmp, false);                
                chi2 = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, stepsize, ar.ppl.xFit_tmp, xstd);               
                if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                    fprintf('Check the Profile at t=%d for inconsistencies, diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_95 + 0.5));
                end
            end
        end       

        if(doPPL)
            ppl_vpl = 'ppl';
            ar.model(m).(data_cond)(c).ppl.(['x' high_low])(jt+t_dir*1, jx)=xSim;
        else
            ppl_vpl = 'vpl';
            ar.model(m).(data_cond)(c).ppl.(['x' high_low '_' ppl_vpl])(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
        end

        ar.ppl.chi2_tmp = chi2;
        ar.model(m).(data_cond)(c).ppl.([ppl_vpl high_low])(jt+t_dir*1, jx)=chi2;
        ar.model(m).(data_cond)(c).ppl.t(jt+t_dir*1,jx)=t_tmp;
        ar.model(m).(data_cond)(c).ppl.(['ps' high_low])(jt+t_dir*1, jx,:)=ar.p;  

        if(t_dir==-1 && jt>1)
            jt=jt-2;
        end
        if(t_dir==-1 && (jt<=1 || t_tmp<ar.model(m).tLim(1) || t_tmp<ar.ppl.options.tEnd))
            fprintf('Backward integration stopped because lower time bound hit, proceeding with normal integration \n');
            break
        end
        if(t_dir == 1 && t_tmp > ar.ppl.options.tEnd)
           fprintf('Hit upper time threshold \n');
           break
        end
    end
    %write LB/UB in ar struct
    if(dir==1)
        struct_string = 'FineUB';
        ppl_string = 'x_high';
    else
        struct_string = 'FineLB';
        ppl_string = 'x_low';
    end
    if(~doPPL)
        ppl_string = [ppl_string '_vpl'];
    end
    if(takeY)
        data_cond = 'data';
        struct_string = ['y' struct_string];
    else
        data_cond = 'condition';
        struct_string = ['x' struct_string];
    end
    ar.model(m).(data_cond)(c).(struct_string)(ar.model(m).(data_cond)(c).tFine<=max(ar.model(m).(data_cond)(c).ppl.t(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx)),jx) = ...
            interp1(ar.model(m).(data_cond)(c).ppl.t(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
            ar.model(m).(data_cond)(c).ppl.(ppl_string)(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
            ar.model(m).(data_cond)(c).tFine(ar.model(m).(data_cond)(c).tFine<=max(ar.model(m).(data_cond)(c).ppl.t(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx))),...
            'pchip',NaN);    
  
    function [xFit_par] = getxFit(xFit)
        xFit_par = xFit;
        x1=ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(it_orig,jx);            
        xFit_par = xFit + (x1-xFit)/(ar.model(m).(data_cond)(c).tFine(it_orig)-t_tmp)*t_dir*stepsize;     
        if(isnan(xFit_par))
           ar.p=pReset;
           if(takeY)
                arLink(true,ar.model(m).(data_cond)(c).tExp(1),true,jx, c, m,NaN);
           end
            fprintf('ERROR IN STEP AT T=%d \n', t_tmp);
            return;
        end
    end   
end

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');
end