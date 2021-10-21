function [chi2_new,par_new] = mergeprofiles(idpar,ple1,ple2)

% [chi2_new,par_new] = mergeprofiles(idpar,ple1,ple2)
%
% This function merges both old and new profile to find the "true" profile,
% which is the minimum of both profile values for every parameter.

q1_notnan = ~isnan(ple1.chi2);
q2_notnan = ~isnan(ple2.chi2);

pars = {ple1.par(q1_notnan,:),ple2.par(q2_notnan,:)};
ps = {ple1.par(q1_notnan,idpar),ple2.par(q2_notnan,idpar)};
chi2s = {ple1.chi2(q1_notnan),ple2.chi2(q2_notnan)};

[soft_lb,id_ub_profile] = max([min(ps{1}),min(ps{2})]);
[soft_ub,id_lb_profile] = min([max(ps{1}),max(ps{2})]);
% soft_lb/ub is the minimal/maximal parameter value of the profile with the
% larger/smaller minimal/maximal value.

q_down = ps{id_lb_profile} < soft_lb;
q_up = ps{id_ub_profile} > soft_ub;

% Interpolate both profiles onto the same parameter values. 
% Note the exeception that profiles can not be extrapolated based on a
% single data point.
ps_unique = unique(round([ps{1};ps{2}],6,'significant'));
if length(ps{1}) > 1
    chi2s_1 = interp1(ps{1},chi2s{1},ps_unique);
    pars_1 = interp1(ps{1},pars{1},ps_unique);
elseif length(ps{1}) == 1
    chi2s_1 = NaN(length(ps_unique),1);
    pars_1 = NaN(length(ps_unique),size(ple1.par,2));
    [~,ind_tmp] = min(abs(ps{1} - ps_unique));
    chi2s_1(ind_tmp) = chi2s{1};
    pars_1(ind_tmp,:) = pars{1};
end
if length(ps{2}) > 1
    chi2s_2 = interp1(ps{2},chi2s{2},ps_unique);
    pars_2 = interp1(ps{2},pars{2},ps_unique);
elseif length(ps{2}) == 1
    chi2s_2 = NaN(length(ps_unique),1);
    pars_2 = NaN(length(ps_unique),size(ple2.par,2));
    [~,ind_tmp] = min(abs(ps{1} - ps_unique));
    chi2s_2(ind_tmp) = chi2s{2};
    pars_2(ind_tmp,:) = pars{2};
end
% chi2s_1/chi2s_2 will likely contain NaNs, since both profiles cover a
% different range of profile parameter values.

n = length(ps_unique);

chi2_new = NaN(n,1);
par_new = NaN(n,size(ple1.par,2));

% Compare values at each interpolated point. Note that this happens only at
% the parameter range where both profiles overlap, since otherwise both
% inequalities will contain a NaN value and return logical zero.
for ii = 1:length(ps_unique)
    % Add values of lower profile to new profile:
    if (chi2s_1(ii) < chi2s_2(ii)) 
        chi2_new(ii) = chi2s_1(ii);
        par_new(ii,:) = pars_1(ii,:); 
    elseif (chi2s_1(ii) >= chi2s_2(ii)) 
        chi2_new(ii) = chi2s_2(ii);
        par_new(ii,:) = pars_2(ii,:);
    else
        continue
    end    
end

% Add values where profiles do not overlap:
if sum(q_down) ~= 0
    chi2_new = [chi2s{id_lb_profile}(q_down);chi2_new];
    par_new = [pars{id_lb_profile}(q_down,:);par_new];
end
if sum(q_up) ~= 0
    chi2_new = [chi2_new;chi2s{id_ub_profile}(q_up)];
    par_new = [par_new;pars{id_ub_profile}(q_up,:)];
end

q_notnan_again = ~isnan(chi2_new);
chi2_new = chi2_new(q_notnan_again);
par_new = par_new(q_notnan_again,:);

if length(par_new(:,idpar)) ~= length(unique(par_new(:,idpar)))
    [~,ind_unis] = unique(par_new(:,idpar));
    chi2_new = chi2_new(ind_unis);
    par_new = par_new(ind_unis,:);
end

end