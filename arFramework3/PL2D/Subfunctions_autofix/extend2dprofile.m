function extend2dprofile

% extend2dprofile
%
% Extends single parameter profiles which are sampled correctly but neither
% reached the confidence threshold nor the parameter bounds. Unfinished
% profiles are identified automatically.

global ar

ar_old = arDeepCopy(ar);

try
    
    smooth2d;
    idpar = ar.ple2d.general.idpar;
    
    % q specifies prediction and direction in which profiles must be
    % extended.
    q = FindExtendables;
    
    if sum(q,'all') ~= 0
        fprintf('\n Extend %i profiles... \n',sum(q,'all'));
    elseif sum(q,'all') == 0
        fprintf('\n No profiles need to be extended. \n');
        return
    end
    
    % Serves as temporary storage for new profiles.
    dummy.chi2 = ar.ple2d.raw.chi2;
    dummy.plpar = ar.ple2d.raw.plpar;
    dummy.par = ar.ple2d.raw.par;
    
    [m,d,idpred,iddata] = Set2dDataPlaceholder;
    
    % Loop over all profiles which need to be extended:
    for ii = 1:2
        for jj = 1:size(q,2)
            if q(ii,jj)
                % Change data value:
                ar.model(m).data(d).yExp(iddata,idpred) = ar.ple2d.smooth.yq(jj,1);
                
                % Initialize ple at first or last parameter value:
                ps_tmp = ar.ple2d.raw.par{jj};
                ps_tmp = ps_tmp(~isnan(ps_tmp(:,1)),:);
                if ii == 1
                    ar.p = ps_tmp(1,:);
                else
                    ar.p = ps_tmp(end,:);
                end
                chi2 = UsePLEFit(idpar);
                chi2_new = chi2 - ar.ple2d.raw.chi2_min;
                del_new = chi2_new - min(ar.ple2d.smooth.zq(jj,:))...
                    - ar.ple2d.smooth.offset;
                
                % Extend profile starting from first or last profile point:
                SetPartialPLEInit(del_new,idpar);
                left_right_tmp = 1:2;
                ple(idpar, 100, [],[],[],[],[],(left_right_tmp == ii))
                pleSmooth;
                
                % Merge existing and extended profile:
                ple1.chi2 = ar.ple2d.raw.chi2(:,jj);
                ple1.par = ar.ple2d.raw.par{jj};                
                ple2.chi2 = (ar.ple.chi2s{idpar} ...
                    - ar.ple2d.raw.chi2_min - ar.ple2d.smooth.offset)';
                ple2.par = ar.ple.ps{idpar};                
                [chi2,par] = mergeprofiles(idpar,ple1,ple2);
                
                % Store results in dummy:
                dummy = StoreFixedResults(dummy,chi2,par,jj,idpar);
            end
        end
        
    end
catch exception
    fprintf(['ERROR extend2dprofile: Resetting ar struct.',...
        ' Error message: \n %s \n Line: %s \n'] ,...
        exception.message, sprintf('%i, ',exception.stack.line));
end

ar = ar_old;
ar.ple2d.raw.chi2 = dummy.chi2;
ar.ple2d.raw.plpar = dummy.plpar;
ar.ple2d.raw.par = dummy.par;

end

function q = FindExtendables

% Find profiles which neither reach the confidence threshold nor reach the
% parameter bounds.

global ar

% eps needs to be chosen in accordance with bounds2d (there: eps = 10^-1). 
% If a confidence bound is in the epsilon neighborhood of the hard 
% parameter bounds, it is set equivalent to the bound. 
% Thus profiles should be at least extended into this neighborhood such that 
% all profiles have defined boundaries.

eps = 10^-2;
idpar = ar.ple2d.general.idpar;
q = zeros(2,length(ar.ple2d.raw.predsteps));

for ii = 1:length(ar.ple2d.raw.predsteps)
    
    chi2s_tmp = ar.ple2d.raw.chi2(:,ii);
    plpars_tmp = ar.ple2d.raw.plpar(~isnan(chi2s_tmp),ii);
    chi2s_tmp = chi2s_tmp(~isnan(chi2s_tmp));
    
    if (chi2s_tmp(1) - min(chi2s_tmp) < ...
        ar.ple2d.config.autofix.excess_factor * icdf('chi2',0.95,...
            ar.ple2d.config.bounds.alpha_bounds)) && ...
            (plpars_tmp(1) - ar.lb(idpar) > eps)
        q(1,ii) = 1;
    end
    if (chi2s_tmp(end) - min(chi2s_tmp) < ...
            ar.ple2d.config.autofix.excess_factor * icdf('chi2',0.95,...
            ar.ple2d.config.bounds.alpha_bounds)) && ...
            (ar.ub(idpar) - plpars_tmp(end) > eps)
        q(2,ii) = 1;
    end
    
end

end

