% individually releases fixed parameters, refits and checks the reduction
% of chi^2
%
% arTrackEffects

% function arTrackEffects

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

stop = false;
pLabelRelease = {};

while(~stop)
    pReset = ar.p;
    qFitReset = ar.qFit;
    chi2Last = arGetMerit('chi2');
    
    chi2yield = zeros(1,sum(ar.qFit == 0));
    jks = find(ar.qFit == 0);
    
    count = 1;
    for j=jks
        fprintf('checking chi^2 yield for fitted %s...\t', ar.pLabel{j});
        
        ar.qFit(j) = 1;
        arFit(true);
        chi2yield(count) = chi2Last - arGetMerit('chi2');
        
        fprintf('%f\n', chi2yield(count));
        
        count = count + 1;
        
        ar.p = pReset;
        ar.qFit = qFitReset;
    end
    fprintf('\n');
    
    labeltmp = ar.pLabel(jks);
    [chi2yieldsorted, ijks] = sort(chi2yield);
    
    [pval,qval] = arLRT(chi2yieldsorted, 1, ar.config.alpha_level);
    
    for j=1:length(chi2yieldsorted)
        if(~qval(j))
            tmpstr = '*';
        else
            tmpstr = ' ';
        end
        
        fprintf('Dchi^2 = %f\t (pval = %f, alpha = %f, %s) for %s\n', ...
            chi2yieldsorted(j), pval(j), ar.config.alpha_level, tmpstr, labeltmp{ijks(j)});
    end
    
    if(sum(~qval)<1)
        stop = true;
        fprintf('\nNo more significant improving parameters found.\n');
        fprintf('Released parameters:\n');
        disp(pLabelRelease')
    else
        fprintf('\nReleasing parameter %s\n\n', ar.pLabel{jks(ijks(end))});
        pLabelRelease{end+1} = ar.pLabel{jks(ijks(end))}; %#ok<SAGROW>
        ar.qFit(jks(ijks(end))) = 1;
        arFit(true);
    end
end

