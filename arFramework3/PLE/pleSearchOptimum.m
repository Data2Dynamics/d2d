% pleSearchOptimum(jks)
% search for better optima, based on down-jump of the profile likelihood
%
% jks  parameters   [1:length(ar.ple.ps)]

function pleSearchOptimum(jks)

global ar

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~exist('jks','var'))
    jks = 1:length(ar.ple.ps);
end

currentChi2 = ar.ple.merit;
currentP = ar.ple.p;

fprintf('PLE actual best chi^2 = %f\n', currentChi2);

trialP = zeros(0,length(currentP));
sourcePLE = [];

dchi2 = 0.5;
dchi2b = 0.5;

for jk=jks
    if(~isempty(ar.ple.chi2s{jk}) && ar.ple.IDstatus(jk)~=4)
        candidateindex = [];
        
        n = length(ar.ple.chi2s{jk});
        chi2s = ar.ple.chi2s{jk};
        
        % down
        for j=ceil(n/2):-1:2
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(1:j));
            if(minindex ~= j && ~isnan(minchi2) && ...
                    (crumin-minchi2)>dchi2 && abs(currentChi2-minchi2)>dchi2b)
				candidateindex = union(candidateindex, minindex); %R2013a compatible
            end
        end
        
        % up
        for j=ceil(n/2):1:n-1
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(j:n));
            if((minindex+j-1) ~= j && ~isnan(minchi2) && ...
                    (crumin-minchi2)>dchi2 && abs(currentChi2-minchi2)>dchi2b)
				candidateindex = union(candidateindex, (minindex+j-1)); %R2013a compatible
            end
        end        
        trialP = [trialP; ar.ple.ps{jk}(candidateindex,:)];
        sourcePLE = [sourcePLE ones(1,length(candidateindex))*jk];
        
%         disp([jk length(candidateindex)])
%         ar.ple.ps{jk}(candidateindex,jk)'
%         ar.ple.chi2s{jk}(candidateindex)
    end
end

trialChi2s = [];

fprintf('PLE tries %i parameter sets:\n', size(trialP,1));

for j = 1:size(trialP,1)
    fprintf('%i/%i\tPLE of %s at %f...', j, size(trialP,1), ...
        ar.ple.p_labels{sourcePLE(j)}, trialP(j,sourcePLE(j)));
    
    try
        feval(ar.ple.integrate_fkt, trialP(j,:));
        p = feval(ar.ple.fit_fkt);
        trialP(j,:) = p;
        trialChi2s(j) = feval(ar.ple.merit_fkt); %#ok<*AGROW>
        fprintf(' chi^2 = %f\n', trialChi2s(j));
    catch %#ok<*CTCH>
        trialP(j,:) = nan(size(currentP));
        trialChi2s(j) = nan;
        fprintf(' FIT ERROR\n');
    end
end

ar.ple.localOpts = [currentP; trialP];
ar.ple.localOptsChi2s = [currentChi2 trialChi2s];

[dummy,imin] = min(ar.ple.localOptsChi2s);
feval(ar.ple.setoptim_fkt, ar.ple.localOpts(imin,:));
fprintf('PLE selected best parameter values\n');
fprintf('PLE new best chi^2 = %f\n', feval(ar.ple.merit_fkt));
