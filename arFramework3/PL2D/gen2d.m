function gen2d(idpar)
% gen2d(idpar)
%
% Generates the 2d-Profile for the specified parameter and experimental
% condition. Call Init2DPL and calculate the validation profile by
% calling InitVPL and VPL to specify the experimental condition 
% before using this function. Make sure an optimum is reached
% before use of this method.
%
% idpar: index of parameter of 2d-profile (mandatory input)
%
% See also: Init2DPL, autofix2dprofile, score2d, bounds2d


global ar

if ~isfield(ar,'ple2d')
    disp('ERROR gen2d: Intitialize 2d profile generation by calling Init2DPL')
    return
end
if ~exist('idpar','var') || isempty(idpar)
    disp('ERROR gen2d: Pass the index of the parameter of interest to gen2d')
end

% Take information from the validation profile as it has to be calculated
% beforehand anyway:
m = ar.vpl.general.m;
d = ar.vpl.general.d;
idpred = ar.vpl.general.idpred;
tpred = ar.vpl.general.tpred; 
sigma = ar.vpl.general.sigma;

% Add complete information of added data point to ar:
ar.ple2d.general.model = m;
ar.ple2d.general.condition = d;
ar.ple2d.general.idpar = idpar;
ar.ple2d.general.pLabel = ar.pLabel{idpar};
ar.ple2d.general.idpred = idpred;
ar.ple2d.general.yLabel = ar.model(m).data(d).y{idpred};
ar.ple2d.general.tpred = tpred;
ar.ple2d.general.sigma = sigma;

ar_old = arDeepCopy(ar);
% Needed to revert temporary changes to ar after function use
npar  = ar.ple2d.config.gen2d.npar;

try    
    fprintf('\n 2D profile for parameter %s and observable %s ... \n',...
        ar.pLabel{idpar},ar.model(m).data(d).y{idpred})

%% Initialize the starting prediction 

    % Starting data point is current optimal prediction:
    ar.model(m).data(d).tExtra = tpred; 
    % Adds new time point to simulated pointd after calling arLink
    arLink(true); 
    arSimu(true,true); % Integrate ODE for fine time grid with new time point
    opt_pred = ar.model(m).data(d).yFineSimu(...
        ar.model(m).data(d).tFine == tpred,idpred); 
    arAddToData(m,d,idpred,tpred,opt_pred,2,1); 
    arLink(true); % Add Data and Link separately to suppress output
    iddata = size(ar.model(m).data(d).yExp,1); %index for new data point
    ar.model(m).data(d).yExpStd(iddata,idpred) = sigma;
    arFit(true); % Fit again because new data point changes optimum
    arSimu(true,true);
    opt_pred = ar.model(m).data(d).yFineSimu(...
        ar.model(m).data(d).tFine == tpred,idpred); 


%% Find predictions at which parameter profiles are sampled

   [predsteps,levels] = findsamplepoints(ar.vpl.results.z',ar.vpl.results.chi2',...
        ar.ple2d.config.gen2d.levels,sigma,ar.ple2d.config.gen2d.weakness); 
    if max(levels) < icdf('chi2',ar.ple2d.config.bounds.alphavpl_max,1)
        fprintf(['WARNING gen2d: The maximal validation profile value ',... 
            '(combined with the weak prior) is %0.4g and thus it is smaller \n',...
            'than the averaging range specified in ',...
            'ar.ple2d.config.bounds.alpha_max. \n'],max(levels))
    end
    % The optional arguments sigma and weakness add a weak prior to the
    % validation profile, thus even flat validation profiles have a finite
    % sampling range        
    
    % Check whether optimal prediction is in accordance with validation
    % profile:
    [~,ind_tmp] = min(abs(opt_pred-predsteps));
    if abs(levels(ind_tmp)-min(levels)) > 1e-4
        fprintf(['\n ERROR gen2d:',...
            ' Current optimal prediction found does not coincide with',...
            ' validation profile minimum. \n Try to initialize 2d-profile',...
            ' calculation from an optimum. \n']);
        ar = ar_old;
        return
    end
    
    predsteps = sort([predsteps,opt_pred]);
    opt_pred_ind = find(predsteps == opt_pred);
    if opt_pred_ind ~= 1
        levels = [levels(1:(opt_pred_ind-1)),0,levels(opt_pred_ind:end)];
    end
    npred = length(predsteps);
    
%% Calculate the parameter profiles for each prediction

    fprintf('\n Calculate Parameter profiles for %i predictions \n', npred)
 
    chi2 =  NaN(2*npar+1,npred);
    plpar = NaN(2*npar+1,npred);
    par = cell(1,npred);
    q = NaN(1,npred);
    
    indexlist = cell(1,2);
    indexlist{1} = opt_pred_ind:npred;
    indexlist{2} = sort(1:(opt_pred_ind-1),'descend');

    for jj = 1:2
        % Determine whether upwards direction has been completed
        if(exist('ii', 'var'))
            % Reset ar.p to starting optimum
            ar.p = par{indexlist{1}(1)}(npar+1,:);
        end
       
        for ii = indexlist{jj}
            if jj == 1
                fprintf('\n Calculate profile number %i \n',ii-(opt_pred_ind-1))
            else
                fprintf('\n Calculate profile number %i \n',npred+1-ii);
            end
            
            % Initialize validation data point:
            ar.model(m).data(d).yExp(iddata,idpred) = predsteps(ii);
            
            % Fit before initializing ple to start calculation in new optimum. 
            % Since the merit function value only changes slightly when 
            % adding a data point, this should still be the 'same' optimum.
            arFit(true);
            % Calculate the current profile:
            arPLEInit(true,true,ar.ple2d.config.gen2d.plmode,0);
            ar.ple.showCalculation = 0;
            ar.ple.alpha_level = 1-ar.ple2d.config.bounds.alpha_bounds;
            ar.ple.dchi2_point = icdf('chi2',...
                ar.ple2d.config.bounds.alpha_bounds,1);
            try
                ple(idpar,npar,ar.ple2d.config.gen2d.relchi);
                pleSmooth(idpar,1); % should automatically smoothen all jumps
                % Store single profiles:
                q(ii) = true;
            catch    
                fprintf('\n WARNING: Profile is skipped \n');
                q(ii) = false;
            end
            chi2(:,ii)  = ar.ple.chi2s{idpar};
            plpar(:,ii) = ar.ple.ps{idpar}(:,idpar);
            par{ii}     = ar.ple.ps{idpar}; 
            
            try
                
                % Check whether new chi2 minimum for the current prediction
                % coincides with the validation profile one:                
                not_nan_vpl = ~isnan(ar.vpl.results.z);
                par_vpl = interp1(ar.vpl.results.z(not_nan_vpl),...
                    ar.vpl.results.ps(not_nan_vpl,idpar),predsteps(ii));
                eps = (10^-8)*(ar.ub(idpar)-ar.lb(idpar));
                if par_vpl >= ar.ub(idpar) - eps
                    par_vpl = ar.ub(ipdar) - eps;
                elseif par_vpl <= ar.lb(idpar) + eps
                    par_vpl = ar.lb(idpar) + eps;
                end
                % par_vpl is the profile parameter value which optimized 
                % the validation profile. 
                
                chi2s_single = chi2(:,ii) - min(chi2(:,ii));
                not_nan_ple = ~isnan(chi2s_single);
                chi2_valmin = interp1(plpar(not_nan_ple,ii),...
                    chi2s_single(not_nan_ple),par_vpl);
                % Get parameter profile value at par_vpl and compare to
                % parameter profile minimum
                
                % If minimum found by parameter profile differs by at least
                % some threshold value, this is saved:
                if chi2_valmin > 0.1
                    if ~exist('maximal_difference','var')
                        maximal_difference = chi2_valmin;
                    else
                        maximal_difference = max([maximal_difference,chi2_valmin]);
                    end
                end
                
            catch
                continue
            end
        end
    end
    chi2_min = min(min(chi2));
    chi2 = chi2 - chi2_min;
    % Shift global optimum to zero
    
    ar = ar_old;
    
    if exist('maximal_difference','var')
        fprintf(['WARNING gen2d: There have been one or more parameter',...
            ' profiles which have a better optimum than indicated',...
            ' by the independent \n validation profile calculation.',...
            ' The maximal difference between the optima is %0.4g \n',...
            ' Set ar.ple2d.config.bounds.vpl_mode = 2',...
            ' to obtain a more accurate score later. \n'],...
            maximal_difference);
    end
    
catch exception
    fprintf(['ERROR gen2d: Resetting ar struct.',...
        ' Error message: \n %s \n Line: %s \n'],...
        exception.message, sprintf('%i, ',exception.stack.line));
    ar = ar_old;
    return
end

%Add raw data of the 2d profile to ar:
q = logical(q);
ar.ple2d.raw.chi2 = chi2(:,q);
ar.ple2d.raw.plpar = plpar(:,q);
ar.ple2d.raw.par = par(q);
ar.ple2d.raw.predsteps = predsteps(q);
ar.ple2d.raw.levels = levels(q);
ar.ple2d.raw.chi2_min = chi2_min;

end
