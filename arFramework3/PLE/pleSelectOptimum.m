% select best optima, based on the profile likelihood

function pleSelectOptimum(jks)

global ar

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~exist('jks','var'))
    jks = 1:length(ar.ple.ps);
end

chi2best = ar.ple.merit + 0;
pbest = [];

fprintf('PLE actual best fit: %f\n', chi2best);

for jk=jks
    if(~isempty(ar.ple.chi2s{jk}))
        [minchi2, minindex] = min(ar.ple.chi2s{jk});
        
        if(minchi2 < chi2best)
            fprintf('PLE #%3i (%20s) at %+f, better fit: %f (%f)\n', jk, ar.ple.p_labels{jk}, ...
                ar.ple.ps{jk}(minindex,jk), minchi2, minchi2 - chi2best);
            pbest = ar.ple.ps{jk}(minindex,:) + 0;
            chi2best =  minchi2 + 0;
        end
    end
end

if(~isempty(pbest))
    feval(ar.ple.setoptim_fkt, pbest);
    fprintf('PLE selected best parameter values...\n');
    fprintf('PLE new best fit: %f\n', ar.ple.merit);
else
    fprintf('PLE not better fit found.\n');
end