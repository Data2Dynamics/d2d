function [chi2s,zs] = VPL
% [chi2s,zs] = VPL
%
% Calculate the validation profile (or prediction profile) for a specific condition.
% Intialize by calling InitVPL. Turn Bessel correction off before use.
% Saving the ar struct before use is recommended, although VPL should return  
% to the previous ar-version if an error occurs (Exception: termination by
% force).
%
% Setting ar.vpl.config.prediction = true switches the exit condition of VPL
% to calculate the prediction profile up to a certain threshold.
% Calculation of prediction profiles works best if sigma is of about the
% same size as the prediction uncertainty.
%
% Visualize results by calling vplPlot.
%
% See also: InitVPL, vplPlot, vplSmooth

global ar

ar_old = arDeepCopy(ar);

if ~isfield(ar,'vpl')
    disp('Warning: Intitialize validation profile calculation by calling InitVPL')
    return
end
if ar.config.useFitErrorCorrection ~= 0
    fprintf(['Warning: Bessel correction is turned on. \n'...
        'Recalculate global optimum and profile with ar.config.useFitErrorCorrection = 0 \n'])
    return
end

try
%% Set up general struct gen_struct to store algorithm results
% This struct will be filled by this algorithm and saved at the end.

    maxsteps = ar.vpl.config.maxsteps;
    m = ar.vpl.general.m;
    d = ar.vpl.general.d;
    idpred = ar.vpl.general.idpred;
    tpred = ar.vpl.general.tpred;
    sigma = ar.vpl.general.sigma;
    
    % Values in the field perm will be saved later
    gen_struct=struct('perm',[],'temp',[]);
    gen_struct.perm.results.chi2 = NaN(2*maxsteps+1,1);
    gen_struct.perm.results.z = NaN(2*maxsteps+1,1);
    gen_struct.perm.results.pred = NaN(2*maxsteps+1,1);
    gen_struct.perm.results.ppl = NaN(2*maxsteps+1,1);
    gen_struct.perm.results.ps = NaN(2*maxsteps+1,length(ar.p));
    gen_struct.perm.test.chi2limit = NaN(2*maxsteps+1,1);
    gen_struct.perm.test.chi2dif = NaN(2*maxsteps+1,1);
    gen_struct.perm.test.delz = NaN(2*maxsteps+1,1);
 
    % Entries in the field temp will not be saved and are used to pass
    % algorithm-specific information to step choice routines.
    gen_struct.temp.maxstepsize = ar.vpl.config.maxstepsize;
    gen_struct.temp.minstepsize = ar.vpl.config.minstepsize;
    gen_struct.temp.chi2dif_min = ar.vpl.config.chi2dif_min;
    gen_struct.temp.chi2dif_max = ar.vpl.config.chi2dif_max;
    gen_struct.temp.stepfactor = ar.vpl.config.stepfactor;
    
%% Add new data point to the existing ar-struct and calibrate

    % Starting data point is current optimal prediction:
    ar.model(m).data(d).tExtra = tpred; %Adds new time point after calling arLink
    arLink(true); 
    arSimu(true,true); %Integrate ODE for fine time grid with new time point
    z_old = ar.model(m).data(d).yFineSimu(...
        ar.model(m).data(d).tFine == tpred,idpred); 
    arAddToData(m,d,idpred,tpred,z_old,2,1); 
    arLink(true); %Add Data and Link separately to suppress output
    iddata = size(ar.model(m).data(d).yExp,1); %index for new data point
    ar.model(m).data(d).yExpStd(iddata,idpred) = sigma;
    
    % Reoptimize because optimal parameters change after adding data point:
    arFit(true); 
    arCalcMerit(ar.vpl.config.sensi,ar.p(ar.qFit==1)); %second argument needed to suppress output
    gen_struct.temp.iddata = iddata;
    chi2_norm = arGetMerit(true);
    ar.vpl.general.chi2_norm = chi2_norm;
    pred_init = ar.model(m).data(d).yExpSimu(iddata,idpred);
    % Note that new optimal prediction may again be different from added
    % data point
    
    fprintf('\n Start at prediction z = %0.4g \n',z_old)
       
%% Save initial struct values before any step:
    chi2_new = 0; %Minimal objective function value is normalized to zero
    p_init = ar.p; 
    gen_struct.perm.results.z(maxsteps+1) = z_old;
    gen_struct.perm.results.chi2(maxsteps+1) = chi2_new;
    gen_struct.perm.results.ps(maxsteps+1,:) = p_init;
    gen_struct.perm.test.chi2dif(maxsteps+1) = 0;
    gen_struct.perm.test.delz(maxsteps +1) = 0;
    gen_struct.perm.test.chi2limit(maxsteps+1) = 0;
    gen_struct.perm.results.pred(maxsteps+1) = pred_init;
    gen_struct.perm.results.ppl(maxsteps+1) = chi2_new-((z_old-pred_init)/sigma)^2;  
    
%% Iteratively change data point value and fit to optimum:   

    for jj = 1:2 
        %Index corresponds to upward and downward direction
        q_exit = 0; %initialize while loop
        
        %Check whethter upwards direction has already been calculated and make
        %corresponding adjustments
        if jj == 1            
            arWaitbar(0);
            delii = 1; 
            %delii introduces a variable minus sign needed for
            %quantities depending on up or down direction
            gen_struct.temp.delii = delii;
        else
            arWaitbar(-1);
            arWaitbar(0);
            delii = -1;
            gen_struct.temp.delii = delii;
            %Reset parameters into global minimum
            ar.p = p_init;
        end
        
        %Iterate steps:
        ii = maxsteps + 1;
        while (q_exit < ar.vpl.config.chi2max)
            
            ii = ii + delii;
            if delii == 1
                arWaitbar(delii*(ii-maxsteps), ar.vpl.config.maxsteps,...
                    'Calculate validation profile towards upper bound');
            else
                arWaitbar(delii*(ii-maxsteps), ar.vpl.config.maxsteps,...
                    'Calculate validation profile towards lower bound');
            end
            
            %Turn off correction based on the actual chi2 difference for
            %the second step, for which the value of the first step was
            %used, because delz has not been adapted for the first step
            if ii == maxsteps+1+2*delii
                gen_struct.temp.correction_switch = 0;
            else
                gen_struct.temp.correction_switch = 1;
            end
            
            %Propose the step:
            if ii == maxsteps+1+delii 
                %first (unadapted) step from starting point
                p_old = p_init;
                delz = ar.vpl.config.firststep; %default unadapted step
                z_old = gen_struct.perm.results.z(maxsteps+1);
                pred_old = pred_init;
                z_new = z_old+delii*delz;
                chi2_old = 0;   
            else
                %any other step
                z_old = z_new;
                chi2_old = chi2_new;
                pred_old = pred_new;
                gen_struct.temp.pred = pred_old; %scalar value for use in step choice function
                delz = feval(ar.vpl.config.step_fkt,...
                    z_old,chi2dif,delz,sigma,gen_struct); %adaptive step choice
                z_new = z_old + delii*delz;
            end
            
            % Make the data step and fit:
            [chi2_new,chi2dif] = vplFit(chi2_old,z_new,iddata,chi2_norm);
            if chi2dif <0 
                % To avoid uncontrolled steps, switch to a mostly conservative
                % default step size and redo the fit.
                ar.p = p_old;
                delz = ar.vpl.config.firststep;
                z_new = z_old+delii*delz;                
                [chi2_new,chi2dif] = vplFit(chi2_old,z_new,iddata,chi2_norm);
            end 
            
            %chi2limit is theoretical upper limit of step:
            chi2limit = (delii*delz/(sigma^2))*...
                (2*(z_old-pred_old)+delii*delz);               
            pred_new = ar.model(m).data(d).yExpSimu(iddata,idpred);
            p_old = ar.p; %needed to reset the parameters if chi2dif<0
            
            %Save algorithm results in struct
            gen_struct.perm.results.z(ii) = z_new;
            gen_struct.perm.results.chi2(ii) = chi2_new;
            gen_struct.perm.results.ps(ii,:) = ar.p;
            gen_struct.perm.test.chi2limit(ii) = chi2limit;
            gen_struct.perm.test.chi2dif(ii) = chi2dif;
            gen_struct.perm.test.delz(ii) = delz;
            gen_struct.perm.results.pred(ii) = pred_new;
            ppl = chi2_new-((z_new-pred_new)/sigma)^2; %alternative exit condition
            gen_struct.perm.results.ppl(ii) = ppl;
            %Plot ppl over pred to obtain prediction profile. 
            
            if ar.vpl.config.showCalculation == true
                fprintf('\n chi2 = %0.4g, z = %0.4g \n',chi2_new,z_new)
            end
            
            %Stopping criteria:
            if (ii == 2*maxsteps+1) || (ii==1)
                break
            end
            if abs(z_new - gen_struct.perm.results.z(maxsteps+1)) > ar.vpl.config.maxrange
                break
            end
            
            if ar.vpl.config.prediction == true
                q_exit = ppl;
            else 
                q_exit = chi2_new;
            end
        end
        
    end
    arWaitbar(-1);
    
catch exception
    fprintf(['\n ERROR VPL: Resetting ar struct. Temporary results are saved in ar.vpl \n',...
        'Error message: %s \n Line: %s \n'],...
        exception.message, sprintf('%i, ',exception.stack.line));
    ar = ar_old;
    ar.vpl = gen_struct.perm;
    return
end

%% Wrap up results

disp('VPL calculation concluded without error.')

%Add values from 'permanent' struct to ar:
gen_struct.perm.config = ar.vpl.config;
gen_struct.perm.general = ar.vpl.general;
ar = ar_old;
ar.vpl = gen_struct.perm;

zs = ar.vpl.results.z;
chi2s = ar.vpl.results.chi2;

if min(chi2s) < 0
    fprintf('\n WARNING VPL: VPL found a better chi2 value by a decrease of %0.4g \n',...
        -min(chi2s));
end

end

function [chi2_new,chi2dif] = vplFit(chi2_old,z_new,iddata,chi2_norm) 
    %Used to reduce code redundancy
    global ar
    
    ar.model(ar.vpl.general.m).data(ar.vpl.general.d).yExp(iddata,...
        ar.vpl.general.idpred) = z_new;
    arFit(true);
    arCalcMerit(ar.vpl.config.sensi,ar.p(ar.qFit==1));
    chi2_new = arGetMerit(true)-chi2_norm;
    chi2dif = chi2_new - chi2_old;
end