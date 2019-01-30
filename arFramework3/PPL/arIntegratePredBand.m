% arIntegratePredBand(general_struct, dir)
% 
% ppl_general_struct    PPL result struct
% dir               direction
%                   dir=1 => uncertainty in upper direction
%                   otherwise => uncertainty in lower direction
% 
% Written by Helge, tried to be documented by Clemens.

function arIntegratePredBand(general_struct, dir)
    global ar
        
    %Set temporal variables from global and structs   
    t_tmp = general_struct.t(ar.ppl.options.whichT);
    doPPL = ar.ppl.options.doPPL;
    gamma_tmp = ar.ppl.options.gamma;
    stepsize = ar.ppl.options.stepsize;
    ed_steps = ar.ppl.options.ed_steps;
    pReset= general_struct.pReset;
    backward = ar.ppl.options.backward;
    fineInt = ar.ppl.options.fineInt;
    takeY = general_struct.takeY;
    m = general_struct.m;
    c = general_struct.c;
    jx = general_struct.jx;
    
    npre=0;
    ppl_vpl = general_struct.ppl_vpl;
    data_cond = general_struct.data_cond;
    x_y = general_struct.x_y;    
          
    %Set starting values for pars and position
    if(dir == 1)
        high_low = 'upperBand';
    else
        high_low = 'lowerBand';        
    end
    ar.ppl.xFit_tmp = ar.model(m).(data_cond)(c).ppl.band.(['xs_' ppl_vpl high_low])(1, jx);
    ar.p=squeeze(ar.model(m).(data_cond)(c).ppl.band.(['ps_' high_low])(1, jx,:))';
    if(size(ar.p,1)~=1)
        error('what happened with the parameters?')
    end

    if(isnan(ar.ppl.xFit_tmp))
        error('Starting value is NaN, check calculation of profiles! \n')
    end
    t_dir = 1;
    jt = 0;
    if(backward)
        jt = ar.ppl.nsteps-1;
        t_dir = -1;
    end
    
    if(takeY)
        xLabel = myNameTrafo(ar.model(m).data(c).y{jx});        
    else
        xLabel = myNameTrafo(ar.model(m).x{jx});        
    end    
    nsteps=ar.ppl.nsteps;

    tcount = 1;
    arWaitbar(0);
    while jt<nsteps
        jt = jt+1;
        
        [~,it_orig] = min(abs(ar.model(m).(data_cond)(c).tFine-t_tmp-t_dir*stepsize));
        
        if(toc>tcount)                   
            arWaitbar((jt+npre), nsteps, sprintf(['PPL-integration (' high_low ') for %s at t=%g %i/%i'], xLabel, t_tmp, jt, nsteps));
            tcount = tcount + 0.5; % update every half second
        end        

        if(ed_steps==true)              
            %VPL part    
            [dps, dx, gamma_tmp] = arIntStepPPL(t_tmp, general_struct, t_dir*stepsize, gamma_tmp);             
            
            ar.p(ar.qFit==1)=ar.p(ar.qFit==1) + dps'*stepsize;
            ar.ppl.xFit_tmp = ar.ppl.xFit_tmp + dx(1)*stepsize;           
            t_tmp = t_tmp + t_dir*stepsize;              
            
            if(fineInt || mod(jt/(floor(nsteps/10)),1)==0)
                [chi2, xSim] = arPPL_Chi2Corr(general_struct, t_tmp);                
            else
                [chi2, xSim] = arPPL_GetChi2(t_tmp, false, general_struct, t_dir*stepsize,ar.ppl.xFit_tmp);
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
            chi2 = arPPL_intCorrection(general_struct, t_tmp, dir); 

        end

        if((chi2 - ar.ppl.chi2_threshold) > 0.2 || (chi2 - ar.ppl.chi2_threshold) < -0.2)

            i_count=0;
            curr_chi2 = arGetMerit('chi2');
            if((takeY && (curr_chi2-ar.res(ar.ppl.resi_tmp)^2)<ar.ppl.chi2_threshold) || (~takeY && curr_chi2<ar.ppl.chi2_threshold))
                while((chi2 - ar.ppl.chi2_threshold) > 0.2 || (chi2 - ar.ppl.chi2_threshold) < -0.2 )  
                    fprintf('Correction at t=%d \n',t_tmp);
                    chi2 = arPPL_intCorrection(general_struct, t_tmp, dir);

                    i_count = i_count + 1 ;
                    if(i_count>1)
                        break;
                    end
                end
            end
            if((chi2 - ar.ppl.chi2_threshold) > 0.2 || (chi2 - ar.ppl.chi2_threshold) < -0.2 )
                fprintf('Try to restart at t=%d diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_threshold));
                ar.p = pReset;
                arLink(true,0.,takeY,jx, c, m,NaN);
                [ar.ppl.xFit_tmp, ar.p] = arPredictionProfile(t_tmp, general_struct, false, dir, ar.ppl.xFit_tmp);                
                chi2 = arPPL_GetChi2(t_tmp, false, general_struct, t_dir*stepsize,ar.ppl.xFit_tmp);               
                if((chi2 - ar.ppl.chi2_threshold) > 0.2 || (chi2 - ar.ppl.chi2_threshold) < -0.2 )
                    fprintf('Check the Profile at t=%d for inconsistencies, diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_threshold));
                end
            end
        end       

        if(doPPL)
            ar.model(m).(data_cond)(c).ppl.band.(['xs_' ppl_vpl high_low])(jt+t_dir*1, jx)=xSim;
        else
            ar.model(m).(data_cond)(c).ppl.band.(['xs_' ppl_vpl high_low])(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
        end

        ar.model(m).(data_cond)(c).ppl.band.([ppl_vpl 'likelihood_' high_low])(jt+t_dir*1, jx)=chi2;
        ar.model(m).(data_cond)(c).ppl.band.tFine_band(jt+t_dir*1,jx)=t_tmp;
        
        %Comment out this line to store parameters during band calculation (and in PPL_init)
        %ar.model(m).(data_cond)(c).ppl.band.(['ps' high_low])(jt+t_dir*1, jx,:)=ar.p;  

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
    
    %Reset data point
    if(takeY)
        arLink(true,0.,true,jx, c, m,NaN);
    end
    arCalcMerit();
    %write LB/UB in ar struct
    if(dir==1)
        struct_string = 'FineUB';
    else
        struct_string = 'FineLB';
    end
    ppl_string = ['xs_' ppl_vpl high_low];
    
    struct_string = [x_y struct_string];
    
    
    %Fill FineLB/FineUB for plotting of bands
    ar.model(m).(data_cond)(c).(struct_string)(1:length(ar.model(m).(data_cond)(c).tFine),jx) = ...
            interp1(ar.model(m).(data_cond)(c).ppl.band.tFine_band(~isnan(ar.model(m).(data_cond)(c).ppl.band.(ppl_string)(:,jx)),jx),...
            ar.model(m).(data_cond)(c).ppl.band.(ppl_string)(~isnan(ar.model(m).(data_cond)(c).ppl.band.(ppl_string)(:,jx)),jx),...
            ar.model(m).(data_cond)(c).tFine,'pchip',NaN);    
            
  
    function [xFit_par] = getxFit(xFit)
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