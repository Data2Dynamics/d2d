% search for better optima, based on down-jump of the profile likelihood

function pleSearchOptimum(jks)

global pleGlobals;

if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 
if(~exist('jks','var'))
    jks = 1:length(pleGlobals.ps);
end

currentChi2 = pleGlobals.chi2;
currentP = pleGlobals.p;

fprintf('PLE actual best chi^2 = %f\n', currentChi2);

trialP = zeros(0,length(currentP));
sourcePLE = [];

dchi2 = 0.5;
dchi2b = 0.5;

for jk=jks
    if(~isempty(pleGlobals.chi2s{jk}) && pleGlobals.IDstatus(jk)~=4)
        candidateindex = [];
        
        n = length(pleGlobals.chi2s{jk});
        chi2s = pleGlobals.chi2s{jk};
        
        % down
        for j=ceil(n/2):-1:2
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(1:j));
            if(minindex ~= j && ~isnan(minchi2) && ...
                    (crumin-minchi2)>dchi2 && abs(currentChi2-minchi2)>dchi2b)
                candidateindex = union(candidateindex, minindex);
            end
        end
        
        % up
        for j=ceil(n/2):1:n-1
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(j:n));
            if((minindex+j-1) ~= j && ~isnan(minchi2) && ...
                    (crumin-minchi2)>dchi2 && abs(currentChi2-minchi2)>dchi2b)
                candidateindex = union(candidateindex, (minindex+j-1));
            end
        end        
        trialP = [trialP; pleGlobals.ps{jk}(candidateindex,:)];
        sourcePLE = [sourcePLE ones(1,length(candidateindex))*jk];
        
%         disp([jk length(candidateindex)])
%         pleGlobals.ps{jk}(candidateindex,jk)'
%         pleGlobals.chi2s{jk}(candidateindex)
    end
end

trialChi2s = [];

fprintf('PLE tries %i parameter sets:\n', size(trialP,1));

for j = 1:size(trialP,1)
    fprintf('%i/%i\tPLE of %s at %f...', j, size(trialP,1), ...
        pleGlobals.p_labels{sourcePLE(j)}, trialP(j,sourcePLE(j)));
    
    try
        feval(pleGlobals.integrate_fkt, trialP(j,:));
        p = feval(pleGlobals.fit_fkt);
        trialP(j,:) = p;
        trialChi2s(j) = feval(pleGlobals.merit_fkt); %#ok<*AGROW>
        fprintf(' chi^2 = %f\n', trialChi2s(j));
    catch %#ok<*CTCH>
        trialP(j,:) = nan(size(currentP));
        trialChi2s(j) = nan;
        fprintf(' FIT ERROR\n');
    end
end

pleGlobals.localOpts = [currentP; trialP];
pleGlobals.localOptsChi2s = [currentChi2 trialChi2s];

[dummy,imin] = min(pleGlobals.localOptsChi2s);
feval(pleGlobals.setoptim_fkt, pleGlobals.localOpts(imin,:));
fprintf('PLE selected best parameter values\n');
fprintf('PLE new best chi^2 = %f\n', feval(pleGlobals.merit_fkt));
