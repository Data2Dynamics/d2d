% select best optima, based on the profile likelihood

function pleSelectOptimum(jks)

global pleGlobals;

if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 
if(~exist('jks','var'))
    jks = 1:length(pleGlobals.ps);
end

chi2best = pleGlobals.chi2 + 0;
pbest = [];

fprintf('PLE actual best fit: %f\n', chi2best);

for jk=jks
    if(~isempty(pleGlobals.chi2s{jk}))
        [minchi2, minindex] = min(pleGlobals.chi2s{jk});
        
        if(minchi2 < chi2best)
            fprintf('PLE #%3i (%20s) at %+f, better fit: %f (%f)\n', jk, pleGlobals.p_labels{jk}, ...
                pleGlobals.ps{jk}(minindex,jk), minchi2, minchi2 - chi2best);
            pbest = pleGlobals.ps{jk}(minindex,:) + 0;
            chi2best =  minchi2 + 0;
        end
    end
end

if(~isempty(pbest))
    feval(pleGlobals.setoptim_fkt, pbest);
    fprintf('PLE selected best parameter values...\n');
    fprintf('PLE new best fit: %f\n', pleGlobals.chi2);
else
    fprintf('PLE not better fit found.\n');
end