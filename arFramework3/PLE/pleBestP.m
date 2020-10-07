% pleBestP
% 
%   This function can be used to translate a better optimum seen in the
%   profile likelhoods to the model, i.e. to the global variabe 'ar'.
% 
%   It has been written to stop ple.m by hand if a better optimium is seen.
%   After calling this function arFit should be applied to optimize all
%   parameters around the new optimium.
% 
% 
% Example:
% arLoad
% arFit
% arPLEInit
% ple    % after a while, a better optimium occurs -> user press CTRL-C to stop
% pleBestP
% arFit

function res = pleBestP
        
global ar 
res = struct;
if(~isempty(ar.ple) & length(fieldnames(ar.ple))>2)
    best = min([ar.ple.chi2s{:}]);
    for i=1:length(ar.ple.chi2s)
        if(~isempty(ar.ple.chi2s{i}) && ar.qFit(i)==0)
            warning([ar.pLabel{i},': qFit=0, but ple available. It seems that ple-calculation has been stop while calculation the profile for this parameter: qFit is now changed to qFit=1 !'])
            ar.qFit(i) = 1;
        end
        
        indbest = find(ar.ple.chi2s{i}==best(1));
        res.ip_best = i;
        res.ifit_best = indbest;

        if(~isempty(indbest))
%             i
            ar.p = ar.ple.ps{i}(indbest(1),:);
        end
    end
else
    disp('ar.ple is empty')
end
