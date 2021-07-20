%--------------------------------------------------------------------------
% Function that decomposes the model into submodels (either defined by the
% user, or found via optimization) and analyses their identifiability.
%--------------------------------------------------------------------------

function [identifiables,nonidentif,obs_states,unobs_states] = decomp(modelname,opts,model_results_folder,submodels)

    global x f h p u w 
    load(modelname,'x','f','h','ics','known_ics'); % 'p' may have been changed in the main file
    
    %==========================================================================
    % Initialize variables:
    identifiables     = [];  % vector of identifiable parameters.
    nonidentif        = [];  % vector of unidentifiable parameters.
    new_nonid_pars    = [];  % new vector of unidentifiable parameters.
    new_ident_pars    = [];  % new vector of identifiable parameters.
    obs_states        = [];  % observable states.
    unobs_states      = [];  % unobservable states.
    new_obs_states    = [];  % new vector of observable states.
    new_unobs_states  = [];  % new vector of unobservable states.  
    obs_inputs        = [];  % observable inputs.
    unobs_inputs      = [];  % unobservable inputs.
    new_obs_inputs    = [];  % new vector of observable inputs.
    new_unobs_inputs  = [];  % new vector of unobservable inputs.
    
    %==========================================================================
    % States included in the submodels:     
    
    if opts.decomp_user == 1
        % The submodels have been defined by the user.
        numsub = size(submodels,2);
        optim_flag = 0;
    else
        % The submodels will be found by optimization. 
        optim_flag = 1;
        % Optimization settings: 
        numsub             = numel(x);      
        problem.x_U        = ones(1,numsub);
        problem.f          = 'objective_fun';
        optim_opts.maxtime = opts.maxOpttime;    
    end
    
    isub = 1;                            
    while isub < numsub + 1      
        skip_elim = 0;
        if optim_flag == 1
            fprintf('\n\n >>> Finding an optimal submodel including state number %d ...\n',isub);
            problem.x_L         = zeros(1,numsub);
            problem.x_L(isub)   = 1; 
            [Results]           = combin_optim(problem,optim_opts);
            estados             = find(abs(Results.xbest)); 
        else
            estados = submodels{isub};
        end
        numstates = numel(estados);   
        
        %======================================================================
        % States, functions, and outputs of the reduced model:
        xred = sym(zeros(numstates,1));
        fred = sym(zeros(numstates,1));
        for i=1:numstates
            xred(i) = x(estados(i));
            fred(i) = f(estados(i));
        end

        hred = [];
        for i=1:numel(h)
            % find the variables appearing in each output:
            varsh = symvar(h(i)).';
            vars = cell(numel(varsh),1);
            for j=1:numel(varsh)
                vars{j} = char(varsh(j));
            end
            % If any of the selected states is among those variables, include that output:
            continuar = 1;
            j = 1;
            while (continuar == 1) && (j <= numel(xred))
                isstatej = strcmpi(char(xred(j)),vars);
                j = j + 1;
                if  find(isstatej) ~= 0
                    hred = [hred; h(i)];
                    continuar = 0;
                end
            end
        end  

        %==========================================================================
        % Variables (parameters, states, inputs) appearing in the equations of 
        % the selected states and outputs:
        total_expr = [fred; hred];
        variaveis = symvar(total_expr).';                        
        vars = cell(numel(variaveis),1);
        for i=1:numel(variaveis)
            vars{i} = char(variaveis(i));
        end

        % which parameters are among those variables?
        parslist = [];
        for i=1:numel(p)
            ispari = strcmpi(char(p(i)),vars);
            if  find(ispari) ~= 0
                parslist = [parslist p(i)];
            end

        end
        
        % which known inputs are among those variables?     
        ured = [];                                    
        for i=1:numel(u)                              
            is_input_i = strcmpi(char(u(i)),vars);    
            if  find(is_input_i) ~= 0                
                ured = [ured u(i)];                   
            end                                       
        end        
        
        % which unknown inputs are among those variables?     
        wred = [];                                    
        for i=1:numel(w)                              
            is_w_i = strcmpi(char(w(i)),vars);    
            if  find(is_w_i) ~= 0                
                wred = [ured w(i)];                   
            end                                       
        end

        % which states (other than those in xred) are among those variables?
        stateslist = [];            
        otherx = setdiff(x,xred); % elements of x that are not in xred
        for i=1:numel(otherx)
            ispari = strcmpi(char(otherx(i)),vars);
            isouti = 0;
            for ii=1:numel(h) % which of these states are measured?:
                isouti = isouti + strcmpi(char(otherx(i)),char(h(ii)));
            end
            if  find(ispari) ~= 0 % if the variable is an extra state...
                if isouti == 0 % ... AND the state is NOT an output...
                    % (since measured states are known, need not be
                    % considered as extra *unknown* parameters):
                    stateslist = [stateslist otherx(i)];
                end 
            end
        end        
        
        if opts.unidentif == 1
            % The 'parameters' in the submodel are the parameters;
            % the states are considered as inputs:
            pred = transpose(parslist);
        else
            % The 'parameters' in the submodel include the parameters, AND
            % the UNMEASURED states considered as additional parameters:            
            pred = transpose([parslist,stateslist]);            
        end

        %==========================================================================
        % Dimensions of the reduced problem: 
        m       = numel(hred);                  % number of outputs
        n       = numel(xred);                  % number of states
        q       = numel(pred);                  % number of unknown parameters
        nw      = numel(wred);                  % number of unknown inputs
        qreal   = numel(parslist);              % number of unknown parameters (without taking into account the states considered as parameters)
        r       = n+q+nw;                       % number of states + parameters + unknown inputs
             
        %==========================================================================
        % The submodel must contain at least one output:
        if m > 0
        
            %==========================================================================
            % Find which states are unmeasured:     
            saidas = cell(m,1);
            for i=1:m
                saidas{i} = char(hred(i));
            end
            ismeasured = zeros(1,n);
            for i=1:n
                ismeasured(i) = sum(strcmpi(char(xred(i)),saidas));
            end
            unmeas_xred_indices = find((1-ismeasured)); % indices (on xred) of the unmeasured states

            %==========================================================================
            % Which initial conditions are unknown:
            xknown = x(find(known_ics)); %#ok<FNDSB>
            knownics = cell(numel(xknown),1);
            for i=1:numel(knownics)
                knownics{i} = char(xknown(i));
            end
            isknown = zeros(1,n);
            for i=1:n
                isknown(i) = sum(strcmpi(char(xred(i)),knownics));
            end
%             unknown_xred_ics_ind = find(1-isknown); % indices (on xred) of the states with unknown initial conditions

            %==========================================================================
            fprintf('\n\n >>> Analysing identifiability of a submodel containing:\n %d states:\n %s',n,char(xred));
            fprintf('\n %d outputs:\n %s',m,char(hred));
            fprintf('\n %d known inputs:\n %s',numel(ured),char(ured));
            fprintf('\n %d unknown inputs:\n %s',nw,char(wred));
            fprintf('\n %d parameters:\n %s',q,char(pred));

            %==========================================================================
            % Build Identifiability-Observability matrix, On(x)
            lastrank = NaN;
            nd = ceil(r/m-1); % minimum number of Lie derivatives
            increaseLie = 1;
            fprintf('\n\n >>> Building the observability-identifiability matrix requires at least %d Lie derivatives',nd);
            fprintf('\n     Calculating derivatives: ');
            tic
            
            %======================================================================
            % Input derivatives:
            clear input_der
            %- Create array of known inputs and set certain derivatives to zero:
            if numel(ured)>0
                for ind_u=1:numel(ured) % create array of derivatives of the inputs
                    input_der(ind_u,:) = [ured(ind_u),sym(strcat(char(ured(ind_u)),sprintf('_d')),[1 nd-1])];  %[1 max(opts.nnzDerU)]
                    input_der(ind_u,(opts.nnzDerU(ind_u)+2):end)=0;
                end 
            else
                input_der = [];
            end
            syms zero_input_der_dummy_name
            
            %- Create array of unknown inputs and set certain derivatives to zero:
            if numel(wred)>0 
                clear w_der
                for ind_w=1:numel(wred) 
                    w_der(ind_w,:) = [wred(ind_w),sym(strcat(char(wred(ind_w)),sprintf('_d')),[1 nd+1])];  % [1 max(opts.nnzDerW)]
                    w_der(ind_w,(opts.nnzDerW(ind_w)+2):end)=0;
                end 
                wlvector_red = reshape(w_der(:,1:end-1),[],1);  % reshape array as a column vector
                wlvector_dot = reshape(w_der(:,2:end),[],1);    % vector of derivatives of the unknown inputs 
                %-- Include as states only nonzero inputs/derivatives:
                [nzi,nzj,nz_wlvec] = find(wlvector_red); %#ok<ASGLU>
                wlvector_red = nz_wlvec;   
                wlvector_dot = wlvector_dot(nzi); 
            else        
                wlvector_red = [];
                wlvector_dot = [];
            end
            
            %======================================================================
            % Augment state vector, dynamics:            
            xaugred = [xred;pred;wlvector_red];                 
            faugred = [fred;zeros(numel(pred),1);wlvector_dot]; 
            
            %======================================================================
            % Define numerical equivalents of the symbolic variables:
            if opts.numeric == 1              
                allvariables = symvar([xaugred;hred;faugred;reshape(input_der,[],1)]);
                numeros = vpa(0.1+rand(size(allvariables)));
            end   
            
            %======================================================================
            % Build Oi:
            onx        = zeros(m*(1+nd),n+q+numel(wlvector_red));
            onx        = sym(onx);
            onx(1:m,:) = jacobian(hred,xaugred);% 1st block
            totaltime  = toc;
            ind        = 0; % Lie derivative index (k)
            lasttime   = 0;                     
            past_Lie   = h;   
            %----------------------------------------------------------------------    
            while ind < nd && lasttime < opts.maxLietime % 2nd and subsequent blocks
                tic
                Lieh = onx((ind*m+1):(ind+1)*m,:)*faugred;                             
                extra_term = 0;
                if ind>0 
                    if numel(ured) > 0
                        for i=1:ind 
                            if i < size(input_der,2) 
                                lo_u_der   = input_der(:,i);   
                                hi_u_der   = input_der(:,i+1);
                                lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                                extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
                            end                        
                        end
                    end
                end       
                past_Lie = Lieh + extra_term;         
                onx(((ind+1)*m+1):((ind+2)*m),:) = jacobian(past_Lie,xaugred); 
                lasttime = toc;
                totaltime = totaltime + lasttime;
                ind = ind+1;
                fprintf('%d ',ind); % index of the last derivative to be computed
            end
            if ind == nd 
                while increaseLie == 1
                    fprintf('\n >>> Observability-Identifiability matrix calculated with %d Lie derivatives',nd);
                    fprintf('\n     (calculated in %d seconds)',totaltime);

                    %==========================================================================
                    % Check identifiability by calculating rank:
                    fprintf('\n >>> Calculating rank...');
                    tic
                    % (Optional) Replace known initial conditions:
                    if opts.replaceICs == 1
                        xind = find(known_ics);
                        if size(ics) ~= size(x)
                            ics = transpose(ics);
                        end
                        onx = subs(onx,x(xind),ics(xind));
                    end
                    % (Optional) Convert from symbolic to numeric:
                    if opts.numeric == 1
                        numonx = subs(onx,allvariables,numeros); 
                    else
                        numonx = onx;
                    end
                    rango = double(rank(numonx));
                    fprintf('\n >>> Rank = %d (calculated in %d seconds)',rango,toc);                 
                    if rango == numel(xaugred)%r  % The submodel is full rank 
                        if opts.unidentif == 0
                            obs_states     = xred(unmeas_xred_indices); 
                            new_ident_pars = pred;
                            new_obs_inputs = wred;  
                        end
                        increaseLie = 0;
                    else % With that number of Lie derivatives the array is not full rank.
                        %----------------------------------------------------------
                        % If there are unknown inputs, we may want to check id/obs of (x,p,w) and not of dw/dt:
                        if numel(wred)>0                           
                            [new_ident_pars,new_nonid_pars,new_obs_states,new_unobs_states,new_obs_inputs,new_unobs_inputs] = ...
                                elim_and_recalc(unmeas_xred_indices,rango,numonx,opts,qreal,pred,xred,identifiables,obs_states,obs_inputs,wred);
                            if ( numel(new_ident_pars)==numel(pred) && numel(new_obs_states)==numel(xred) && numel(intersect(wred,new_obs_inputs))==numel(wred) )
                                obs_states    = xred;
                                obs_inputs    = wred;
                                identifiables = pred;
                                increaseLie   = 0;  
                            end       
                        end
                        %----------------------------------------------------------  
                        % If possible (& necessary), calculate one more Lie derivative and retry:
                        if nd < numel(xaugred) && lasttime < opts.maxLietime && rango ~= lastrank && increaseLie == 1
                            tic
                            nd = nd+1;  
                            ind = nd;  
                            extra_term = 0; % reset for each new Lie derivative

                            %-  Known input derivatives: --------------------------
                            if numel(ured) > 0 % Extra terms of extended Lie derivatives
                                % may have to add extra input derivatives (note that 'nd' has grown):
                                clear input_der
                                for ind_u=1:numel(ured) 
                                    input_der(ind_u,:) = [ured(ind_u),sym(strcat(char(ured(ind_u)),sprintf('_d')),[1 nd-1])];
                                    input_der(ind_u,(opts.nnzDerU(ind_u)+2):end)=0;
                                    if opts.numeric == 1 % numerical replacement of the new variable:
                                        if nd-1<=opts.nnzDerU(ind_u)
                                            allvariables = [allvariables,input_der(end)];
                                            numeros      = [numeros,vpa(0.1+rand)];
                                        end
                                    end
                                end                        
                                for i=1:ind-1 
                                    if i < size(input_der,2) 
                                        lo_u_der   = input_der(:,i);   
                                        hi_u_der   = input_der(:,i+1);
                                        lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                                        extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
                                    end 
                                end
                            end
                            
                            %- Unknown input derivatives: ------------------------- 
                            % (add new derivatives, if they are not zero:
                            if numel(wred)>0
                                prev_size = numel(wlvector_red);
                                clear w_der
                                for ind_w=1:numel(wred) 
                                    w_der(ind_w,:) = [wred(ind_w),sym(strcat(char(wred(ind_w)),sprintf('_d')),[1 nd+1])];  
                                    w_der(ind_w,(opts.nnzDerW(ind_w)+2):end)=0;
                                    if opts.numeric == 1 % numerical replacement of the new variable:
                                        if nd+1<=opts.nnzDerW(ind_w)
                                            allvariables = [allvariables,w_der(end)];
                                            numeros      = [numeros,vpa(0.1+rand)];
                                        end
                                    end
                                end                         
                                wlvector_red = reshape(w_der(:,1:end-1),[],1);  % reshape array as a column vector
                                wlvector_dot = reshape(w_der(:,2:end),[],1);    % vector of derivatives of the unknown inputs 
                                %-- Include as states only nonzero inputs/derivatives:
                                [nzi,nzj,nz_wlvec] = find(wlvector_red); %#ok<ASGLU>
                                wlvector_red = nz_wlvec;   
                                wlvector_dot = wlvector_dot(nzi); 
                                %-- Augment state vector & dynamics with new input derivs:
                                xaugred = [xred;pred;wlvector_red];                      
                                faugred = [fred;zeros(numel(pred),1);wlvector_dot]; 
                                %-- Augment size of the Obs-Id matrix if needed:
                                new_size = numel(wlvector_red);
                                onx = [onx,zeros(ind*m,new_size-prev_size)];
                            end  
                            
                            newLie   = onx(((ind-1)*m+1):ind*m,:)*faugred; 
                            past_Lie = newLie + extra_term;        
                            newOnx   = jacobian(past_Lie,xaugred); 
                            onx      = [onx; newOnx];   
                            clear newLie newOnx
                            lasttime  = toc;
                            totaltime = totaltime + lasttime;
                            obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,nd);
                            obsidentfile = strcat(pwd,filesep,'results',filesep,obsidentmatrix);
                            save(obsidentfile);
                            lastrank = rango;  

                        % If that is not possible, there are several possible causes:
                        else
                            if nd >= numel(xaugred) % The maximum number of Lie derivatives has been reached
                                fprintf('\n    The submodel is structurally unidentifiable');
                            else
                                if rango == lastrank
                                    numonx = numonx(1:(end-m),:);
                                    nd = nd - 1;
                                else
                                    if lasttime >= opts.maxLietime
                                        if opts.unidentif == 1
                                            fprintf('\n    => More Lie derivatives would be needed to see if the submodel is structurally unidentifiable.');
                                            skip_elim = 1;
                                            increaseLie = 0;
                                        else
                                            fprintf('\n    => More Lie derivatives would be needed to see if the submodel is structurally identifiable.');
                                        end
                                        fprintf('\n    However, the maximum computation time allowed for calculating each of them has been reached.');
                                        fprintf('\n    You can increase it by changing <<opts.maxLietime>> (currently opts.maxLietime = %d)',opts.maxLietime);

                                    end
                                end
                            end
                            if skip_elim == 0                               
                                % Eliminate columns one by one to check identifiability of the associated parameters:
                                [new_ident_pars,new_nonid_pars,new_obs_states,new_unobs_states,new_obs_inputs,new_unobs_inputs] = ...
                                    elim_and_recalc(unmeas_xred_indices,rango,numonx,opts,qreal,pred,xred,identifiables,obs_states,obs_inputs,wred);
                                increaseLie = 0;
                            end
                        end
                    end
                end

            else% If the maxLietime has been reached, but the minimum of Lie derivatives has not been calculated:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if opts.unidentif == 0
                    fprintf('\n    More Lie derivatives would be needed to analyse the submodel.');
                    fprintf('\n    However, the maximum computation time allowed for calculating each of them has been reached.');
                    fprintf('\n    You can increase it by changing <<opts.maxLietime>> (currently opts.maxLietime = %d)',opts.maxLietime);
                    fprintf('\n >>> Calculating rank...');
                    tic
                    if opts.replaceICs == 1
                        xind = find(known_ics);
                        if size(ics) ~= size(x)
                            ics = transpose(ics);
                        end
                        onx = subs(onx,x(xind),ics(xind));
                    end  
                    if opts.numeric == 1
                        allvariables = symvar(onx);
                        numeros = vpa(0.1+rand(size(allvariables)));
                        numonx = subs(onx,allvariables,numeros); 
                    else
                        numonx = onx;
                    end
                    rango = double(rank(numonx));
                    fprintf('\n  Rank = %d (calculated in %d seconds)',rango,toc); 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [new_ident_pars,new_nonid_pars,new_obs_states,new_unobs_states,new_obs_inputs,new_unobs_inputs] = ...
                        elim_and_recalc(unmeas_xred_indices,rango,numonx,opts,qreal,pred,xred,identifiables,obs_states,obs_inputs,wred); 
                end
            end          
        
        %==========================================================================
        % Add the newly characterized parameters to the vectors of identifiable / non identifiable parameters:
        if size(identifiables,2) ~= 1
            identifiables = transpose(identifiables);
        end
        if size(new_ident_pars,2) ~= 1
            new_ident_pars = transpose(new_ident_pars);
        end              
        identifiables = [identifiables; new_ident_pars];
        identifiables = symvar(identifiables);    

        if size(nonidentif,2) ~= 1
            nonidentif = transpose(nonidentif);
        end
        if size(new_nonid_pars,2) ~= 1
            new_nonid_pars = transpose(new_nonid_pars);
        end
        nonidentif = [nonidentif; new_nonid_pars];
        nonidentif = symvar(nonidentif);

       % And of states / initial conditions:
        if size(new_obs_states,2) ~= 1
            new_obs_states = transpose(new_obs_states);
        end  
        obs_states   = [obs_states; new_obs_states];
        obs_states   = symvar(obs_states);
        if size(new_unobs_states,2) ~= 1
            new_unobs_states = transpose(new_unobs_states);
        end
        unobs_states = [unobs_states; new_unobs_states];
        unobs_states = symvar(unobs_states);
        
        % And of inputs:
        if size(new_obs_inputs,2) ~= 1
            new_obs_inputs = transpose(new_obs_inputs);
        end
        obs_inputs   = [obs_inputs; new_obs_inputs];
        obs_inputs   = symvar(obs_inputs);
        if size(new_unobs_inputs,2) ~= 1
            new_unobs_inputs = transpose(new_unobs_inputs);
        end
        unobs_inputs = [unobs_inputs; new_unobs_inputs];
        unobs_inputs = symvar(unobs_inputs);
                
        end
        
        %==========================================================================
        % Save results:
        provresultsname = sprintf('partial_id_results_%s_submodel_%d',modelname,isub);
        fullprovresultsname = strcat(pwd,filesep,model_results_folder,filesep,provresultsname);
        save(fullprovresultsname);       
        clear Results estados numstates fred hred xred pred wred ured variaveis vars...
            m n q nw r qreal parslist xaugred faugred nd Lieh onx ident rango rangototal
        
        %==========================================================================
        % Proceed with the next submodel:
        isub = isub + 1;  
    end  
end
