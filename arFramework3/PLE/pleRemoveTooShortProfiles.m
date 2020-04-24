% ple = pleRemoveTooShortProfiles(ple,lb,ub)
% 
%   all profiles that are not finished, i.e.  crossed the threshold or
%   reached the bound are eliminated (chi2s is set to NaN) 
% 
%   Example: 
% ple= pleRemoveTooShortProfiles(ple,ar.lb,ar.ub);


function ple = pleRemoveTooShortProfiles(ple,lb,ub)


bool_raus = false(size(ple.conf_ub_point));
bool_raus = bool_raus | isnan(ple.conf_ub_point) | isnan(ple.conf_lb_point);


for i=1:length(ple.ps)
    delta = (ub(i)-lb(i))/20;
    
    if isinf(ple.conf_ub_point(i)) && nanmax(ple.ps{i}(:,i))<(ub(i)-delta)
        bool_raus(i) = true;
    end
    if isinf(ple.conf_lb_point(i)) && nanmin(ple.ps{i}(:,i))>(lb(i)+delta)
        bool_raus(i) = true;
    end
    
end

raus = find(bool_raus);

for j=1:length(raus)
    if length(ple.chi2s)>=raus(j)
        if ~isempty(ple.chi2s{raus(j)})
            ple.chi2s{raus(j)} = NaN*ple.chi2s{raus(j)};
        end
    end
%     ple.ps{raus(j)} = cell(0);
end
