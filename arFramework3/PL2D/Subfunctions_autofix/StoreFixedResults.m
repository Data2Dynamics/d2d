function dummy = StoreFixedResults(dummy,chi2,par,indy,idpar)

% dummy = StoreFixedResults(dummy,chi2,par,indy,idpar)
%
% Increase struct dimensions of ar.ple2d if new profile has more points
% than maximally anticipated.

dim_orig = size(dummy.chi2); % current size
dim_alt = length(chi2); % new size

if dim_orig(1) < dim_alt
    dummy.chi2 = [dummy.chi2;NaN(dim_alt-dim_orig(1),dim_orig(2))];
    dummy.plpar = [dummy.plpar;NaN(dim_alt-dim_orig(1),dim_orig(2))];
    for jj = 1:dim_orig(2)
        dummy.par{jj} = [dummy.par{jj};...
            NaN(dim_alt-dim_orig(1),size(par,2))];
    end
end

% Add new profile to dummy struct:
dummy.chi2(:,indy) = NaN(size(dummy.chi2(:,indy)));
dummy.chi2(1:dim_alt,indy) = chi2;
dummy.plpar(:,indy) = NaN(size(dummy.plpar(:,indy)));
dummy.plpar(1:dim_alt,indy) = par(:,idpar);
dummy.par{indy} = NaN(size(dummy.par{indy}));
dummy.par{indy}(1:dim_alt,:) = par;
        
end

