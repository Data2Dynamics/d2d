% better = pleSelectOptimum(jks)
% 
% select best optima, based on the profile likelihood
% 
%   jks     parameters   [1:length(ar.ple.ps)]
%   better  empty, if no better fit found
%           otherwise struct with fields
% 
%       better.chi2old  chi2 from the fit, i.e. before the profiles are examined
%       better.ip:      which profile has a better chi2
%       better.ifit:    the fit index which shows most improvement
%       better.chi2better: chi2s which are better 
%       better.atEnd:   -1 = left corner, 1 = right corner, 0 = somewhere in the middle
% 
% Example:
% arLoad
% arPLEInit
% ple
% while(~isempty(pleSelectOptimum))
%     pleBestP
%     arFit
%     arPLEInit
%     ple
% end


function better = pleSelectOptimum(jks)

global ar

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~exist('jks','var'))
    jks = 1:length(ar.ple.ps);
end

chi2best = ar.ple.merit + 0;
chi2best0 = chi2best; % best fit without improvement from profiles
pbest = [];

fprintf('PLE actual best fit: %f\n', chi2best);

better = struct;
better.chi2old = chi2best0;
better.ip = [];
better.ifit = [];
better.chi2better = [];
better.atEnd = []; % -1 = left corner, 1 = right corner, 0 = somewhere in the middle

for jk=jks
    if(~isempty(ar.ple.chi2s{jk}))
        [minchi2, minindex] = min(ar.ple.chi2s{jk});
        
        if(minchi2 < chi2best0) % better than the best fit without considering the profiels
            fprintf('PLE #%3i (%20s) at %+f, better fit: %f (%f)\n', jk, ar.ple.p_labels{jk}, ...
                ar.ple.ps{jk}(minindex,jk), minchi2, minchi2 - chi2best0);
            
            better.ip = [better.ip,jk];
            better.ifit = [better.ifit, minindex];
            better.chi2better = [better.chi2better,minchi2];
            if minindex==1
                better.atEnd = [better.atEnd,-1];
            elseif minindex == length(ar.ple.chi2s{jk})
                better.atEnd = [better.atEnd,1];
            else
                better.atEnd = [better.atEnd,0];
            end
            
            if minchi2 < chi2best % better than the previous best
                pbest = ar.ple.ps{jk}(minindex,:) + 0;
                chi2best =  minchi2 + 0;
            end
            

        end
    end
end

if(~isempty(pbest))    
    feval(ar.ple.setoptim_fkt, pbest);
    fprintf('PLE selected best parameter values...\n');
    fprintf('PLE new best fit: %f\n', ar.ple.merit);
    
    [better.chi2best,ibest] = min(better.chi2better);
    better.ip_best = better.ip(ibest);
else
    fprintf('PLE not better fit found.\n');

    better = [];  % make better empty (otherwise it is a struct)
end