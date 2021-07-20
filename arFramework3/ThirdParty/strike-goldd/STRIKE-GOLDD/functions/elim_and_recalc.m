%--------------------------------------------------------------------------
% Function that determines identifiability of individual parameters one 
% by one, by successive elimination of its column in the identifiability 
% matrix and recalculation of its rank
%--------------------------------------------------------------------------
function [new_ident_pars,new_nonid_pars,new_obs_states,new_unobs_states,new_obs_in,new_unobs_in] = elim_and_recalc(unmeas_xred_indices,rangoinicial,numonx,opts,varargin)

global p x unidflag wlvector

switch nargin 
    case 4 
        % called when there is no decomposition nor 'w'      
        pred          = p;
        xred          = x; 
        wred          = wlvector;
        identifiables = [];
        obs_states    = [];
        obs_inputs    = [];
        q             = numel(pred);
        n             = numel(xred);
        nw            = numel(wred);
        qreal         = q;
    case 7 
        % called when there are 'w', but no decomposition
        pred          = p;
        xred          = x; 
        wred          = wlvector;
        identifiables = varargin{1};
        obs_states    = varargin{2};    
        obs_inputs    = varargin{3};
        q             = numel(pred);
        n             = numel(xred);
        nw            = numel(wred); 
        qreal         = q;
    otherwise
        % called in case of decomposition
        % varargin: qreal,pred,xred,identifiables,obs_states,obs_inputs
        qreal         = varargin{1};
        pred          = varargin{2};
        xred          = varargin{3};
        identifiables = varargin{4};
        obs_states    = varargin{5};    
        obs_inputs    = varargin{6};
        wred          = varargin{7};
        q             = numel(pred);
        n             = numel(xred);
        nw            = numel(wred);
end
r  = size(numonx,2); % before: q+n+nw; but with unknown inputs there may also be derivatives
new_ident_pars   = identifiables;
new_nonid_pars   = [];
new_obs_states   = obs_states;
new_unobs_states = [];
new_obs_in       = obs_inputs;
new_unobs_in     = [];
    

%==========================================================================
% ELIMINATE A PARAMETER:
%==========================================================================
% At each iteration we try removing a different column from onx:
for ind=1:qreal % only the first 'qreal' elements of pred are parameters; the following 'q-qreal' are states considered as parameters 
    isidentifiable = ismember(pred(ind),identifiables);
    if isidentifiable
        fprintf('\n Parameter %s has already been classified as identifiable.',char(pred(ind)))
    else 
        indices = 1:r;
        indices(n+ind) = []; % remove the column whose identifiability we want to check
        num_rank = rank(numonx(:,indices));
        if num_rank == rangoinicial
            if (opts.unidentif == 1) || (opts.forcedecomp == 0 && opts.decomp == 0 && unidflag == 1) %%%
            	fprintf('\n    => Parameter %s is structurally unidentifiable',char(pred(ind)));
            	new_nonid_pars = [new_nonid_pars; pred(ind)];
                unidflag = 1;
            else
            	fprintf('\n    => We cannot decide about parameter %s at the moment',char(pred(ind)));
            end  
        else
            if opts.unidentif == 0
            	fprintf('\n    => Parameter %s is structurally identifiable',char(pred(ind)));
            	new_ident_pars = [new_ident_pars; pred(ind)];  
            else
            	fprintf('\n    => We cannot decide about parameter %s at the moment',char(pred(ind)));
            end      
        end
    end
end

if nargin == 4 || 7 % if there is no decomposition
    %==========================================================================
    % ELIMINATE A STATE:
    %==========================================================================
    % At each iteration we try removing a different state from 'xred':
    if opts.checkObser == 1
        for ind=1:numel(unmeas_xred_indices) % for each unmeasured state
            original_index = unmeas_xred_indices(ind); % in this script, 'original_index' refers to xred
            isobservable = ismember(xred(original_index),obs_states);
            if isobservable
                fprintf('\n State %s has already been classified as observable.',char(xred(original_index)))
            else              
                indices = 1:r;
                indices(original_index) = []; %indices(original_index==unmeas_xred_indices) = [];
                num_rank = rank(numonx(:,indices));
                if num_rank == rangoinicial
                    if (opts.unidentif == 1) || (opts.forcedecomp == 0 && opts.decomp == 0 && unidflag == 1) %%%
                        fprintf('\n    => State %s is unobservable',char(xred(original_index)));
                        new_unobs_states = [new_unobs_states; xred(original_index)];
                    else
                        fprintf('\n    => We cannot decide about state %s at the moment',char(xred(original_index)));
                    end 
                else
                    if opts.unidentif == 0
                        fprintf('\n    => State %s is observable',char(xred(original_index)));
                        new_obs_states = [new_obs_states; xred(original_index)];  
                    else
                        fprintf('\n    => We cannot decide about state %s at the moment',char(xred(original_index)));
                    end   
                end
            end
        end
    end

    %==========================================================================
    % ELIMINATE AN UNKNOWN INPUT:
    %==========================================================================
    % At each iteration we try removing a different column from onx:
    for ind=1:nw 
        isobservable = ismember(wred(ind),obs_inputs);
        if isobservable
            fprintf('\n Input %s has already been classified as observable.',char(wred(ind)))
        else 
            indices = 1:r;
            indices(n+q+ind) = []; % remove the column whose identifiability we want to check
            num_rank = rank(numonx(:,indices));
            if num_rank == rangoinicial
                if (opts.unidentif == 1) || (opts.forcedecomp == 0 && opts.decomp == 0 && unidflag == 1) %%%
                    fprintf('\n    => Input %s is unobservable',char(wred(ind)));
                    new_unobs_in = [new_unobs_in; wred(ind)];
                else
                    fprintf('\n    => We cannot decide about input %s at the moment',char(wred(ind)));
                end  
            else
                if opts.unidentif == 0
                    fprintf('\n    => Input %s is observable',char(wred(ind)));
                    new_obs_in = [new_obs_in; wred(ind)];  
                else
                    fprintf('\n    => We cannot decide about input %s at the moment',char(wred(ind)));
                end      
            end
        end
    end
    
end

end
