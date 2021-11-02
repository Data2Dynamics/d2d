function [alpha_pred,pred] = findvalidationpoints
% [alpha_pred,pred] = findvalidationpoints
%
% Calculate the points on the valdiation profile for which boundaray
% averaging is performed. The outputs are the correpsonding predictions and
% alpha levels, and the corresponding parameter profiles are implicitly
% saved in fields in ar.ple2d.bounds. 
%
% Determining the validation alpha levels can pe perfomred from the raw 
% profiles (vpl_mode = 1) or from the smoothed 2d-grid (vpl_mode = 2). 
% The 'only' difference is that by using the second mode, the validation 
% profile is obtained from the 2d-profile instead of the independent 
% calculation. This might be appropriate if the 2d-profile induces a
% different validation profile than the independent validation profile
% calculation.
%
% Inputs for this function are implicitly taken from ar.ple2d.config.bounds.
%
% Outputs:
% pred: Sampled prediction points 
% alpha_pred: Validation confidence levels corresponding to predictions in 
%               pred.
%
% See also: bounds2d, score2d

global ar

nvpl = ar.ple2d.config.bounds.nvpl;

%% Take grid values or unedited profile values
if ar.ple2d.config.bounds.vpl_mode == 1
    par = ar.ple2d.raw.plpar;
    pred = ar.ple2d.raw.predsteps;
    chi2 = ar.ple2d.raw.chi2;
    alphavpls_set = cdf('chi2',ar.ple2d.raw.levels,1);
else
    % For compatibility with the other computation relying only on profiles,
    % parameter and chi2 data must be converted into profile-like data.
    % 1. Find profile values for the levels used in profile calculation by
    %       interpolation
    % 2. Introduce same NaNs
    pred = ar.ple2d.smooth.yq(:,1)';
    chi2 = ar.ple2d.smooth.zq;
    par = ar.ple2d.smooth.xq;
    
    % The grid profiles must first be interpolated onto the levels stored in 
    % ar.ple2d.raw.levels to be compatible.
    level_tmp = unique(round(ar.ple2d.raw.levels,5));
    [~,yvpl,zvpl] = vplfrom2d;
    [pred_new,levels_new] = findsamplepoints(yvpl,zvpl,level_tmp,...
        ar.ple2d.general.sigma,ar.ple2d.config.gen2d.weakness);
    % Find predictions from the 2D-calculated validation profile which
    % correspond to the chi2-levels in ar.ple2d.raw.levels. This imposes the
    % actual validation confidence levels found, not the ones found in the
    % independent validation profile calculation.
    alphavpls_set = cdf('chi2',levels_new,1);
    
    par_new = NaN(size(chi2,2),length(pred_new));
    chi2_new = NaN(size(chi2,2),length(pred_new));
    q_nan = logical(zeros(1,length(pred_new)));
    for ii = 1:length(pred_new)
        k = 1;           
        while (pred_new(ii) > pred(k)) 
            if k == length(pred) 
                break
            else
                k = k+1;
            end
        end     
        if ((k == 1) || (k == length(pred))) && (pred_new(ii) > pred(k))
            % The and-connection can occur if k == length(pred)
            fprintf([' \n WARNING findvalidationpoints: Prediction to certain ',...
                'confidence level yielded an unexpected value. \n Skip this ',...
                'prediction. \n']);
            continue
        else
            %Interpolate profiles:
            q_nan(ii) = true;
            rel_distance = (pred_new(ii)-pred(k-1))/(pred(k)-pred(k-1));
            par_new(:,ii) = par(1,:);
            chi2_new(:,ii) = InterProfiles(par_new(:,ii),...
                chi2(k-1,:),chi2(k,:),rel_distance);
        end
    end    
    pred = pred_new(q_nan);
    par = par_new(:,q_nan);
    chi2 = chi2_new(:,q_nan);
    alphavpls_set = alphavpls_set(q_nan);
        
    chi2_notnan = ~isnan(chi2(:,:));
    chi2_ones = ones(size(chi2,1),size(chi2,2));
    chi2_ones(chi2_notnan == 0) = NaN;
    par = chi2_ones.*par;
    %Introduce the same Nans in par as in chi2
end

%% The core of the function starts here:
% Find the validation data points closest to the values given in alphavpls_guide.

[alphavpls_unique,~,ind_levels_unique] = unique(round(alphavpls_set,5));
alphavpls_guide = linspace(ar.ple2d.config.bounds.alphavpl_min,...
    ar.ple2d.config.bounds.alphavpl_max,nvpl);
alphavpls_guide = [0,alphavpls_guide];
%alphavpls_set: alpha levels of all different parameter profiles
%alphavpls_unique: alphavpls_set has the same values in both directions
%   from the minimum, reduce to unique values
%ind_levels_unique: Indices for the back-transformation
%alphavpls_guide: Prespecified alpha values determining at which levels
%    parameter profiles should be sampled (mainly user input)

%Choose data alpha-levels as close as possible to specified alpha_levels
q_eff = zeros(length(alphavpls_set),1);
for ii = 1:(nvpl+1)
    [~,ind_min_tmp] = min(abs(alphavpls_unique - alphavpls_guide(ii)));
    q_eff = (q_eff | (ind_levels_unique == ind_min_tmp));
end
[ind_levels_unique2,q_rmneighbours] = removeneighbours(ind_levels_unique(q_eff));
alpha_pred = alphavpls_unique(ind_levels_unique2);
%q_eff: Logical vector containing the indices of all predictions used for averaging
%q_rmneighbours: Logical vector removing neighbouring values in the index
%   vector, so that there is only one set of bounds corresponding to each
%   confidence level
%alpha_pred: Resulting confidence levels with the exact same values for
%   repetitions

%Reduce data sets to relevant profiles:
pred = pred(q_eff);
pred = pred(q_rmneighbours);
par = par(:,q_eff);
par = par(:,q_rmneighbours);
chi2 = chi2(:,q_eff);
chi2 = chi2(:,q_rmneighbours);

ar.ple2d.bounds.pred_points = pred;
ar.ple2d.bounds.alpha_pred = alpha_pred;
ar.ple2d.bounds.pars_profile = par;
ar.ple2d.bounds.chi2s_profile = chi2;
end

function [ind_new,q] = removeneighbours(ind)
% Removes neighbouring values in an integer vector. q is the logical vector
% determining which values remain, while ind_new is the index of the
% these values

q = true;
for ii = 1:(length(ind)-1)
    if ind(ii+1)-ind(ii) == 0
        q(ii+1) = false;
    else
        q(ii+1) = true;
    end
end
ind_new = ind(q);

end

function z_new = InterProfiles(x,z1,z2,rel_distance)
% x is profile parameter value
% y is prediction index (rel_distance is the normalized distance) 
% z is profile value

% Bring everything into the right form:
x_data = [x,x];
y_data = [ones(size(x)),2*ones(size(x))];
z_data = [z1,z2];
q = ~isnan(z_data);
x_data = x_data(q);
y_data = y_data(q);
z_data = z_data(q);
xq = x;
yq = (1+rel_distance)*ones(size(x));

% Interpolate between both profiles:
z_new = griddata(x_data,y_data,z_data,xq,yq);


end
