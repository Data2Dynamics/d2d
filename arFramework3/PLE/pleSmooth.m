% pleSmooth(jk, [quick], [point], [dr])
% 
% Smoothen jumps (peaks) in a likelihood profile which might occure because
% of local optima. Removing these peaks is performed by repeating
% optimization by other initial guesses (for the non-profiled parameters).
%
%   jk         parameter index for the profile to be smoothend
%   quick      skip selection and pick the first trial point 
%              [default = 0]
%   point      provide custom parameter value to start at (closest point
%              will be chosen)
%   dr         initial search direction (-1 or 1; required when point is 
%              specified)
% 
% Sometimes it is reasonable to call pleSmooth repeatedly to eliminate all
% peaks.
% 
% Example:
% arPLEInit
% ple(1)
% pleSmooth

function pleSmooth(jk, quick, point, dr)

global ar

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~isfield(ar.ple, 'showCalculation'))
    ar.ple.showCalculation = true;
end
if(~exist('quick','var'))
    quick = false;
end
if(nargin<1)
    fprintf('PLE smoothing for %i parameters ...\n', sum(ar.qFit))
    jindex = find(ar.qFit);
    for j=1:length(ar.ple.chi2s)
        pleSmooth(jindex(j))
    end
    return
elseif(length(jk)>1)
    fprintf('PLE smoothing for %i parameters ...\n', length(jk))
    for j=1:length(jk)
        pleSmooth(jk(j))
    end
    return
end

if (nargin>3)
    if (nargin<4)
        error('Need to specify point and direction');
    end
end

dchi2 = 0.1;

if(~ar.qFit(jk))
    fprintf('\nPLE#%i smoothing SKIPPED: parameter %s is fixed\n', jk, ar.ple.p_labels{jk});
    return;
elseif(jk <= length(ar.ple.chi2s) && ~isempty(ar.ple.chi2s{jk}) && ar.ple.IDstatus(jk)~=4)
    fprintf('\nPLE#%i smoothing for parameter %s\n', jk, ar.ple.p_labels{jk});
    
    
    trialP = zeros(0,length(ar.ple.p));

    n = length(ar.ple.chi2s{jk});
    chi2s = ar.ple.chi2s{jk};
    [~, globminindex] = min(chi2s);
    
    if ( nargin < 3 )
        candidateindex_down = [];
        candidateindex_up = [];
        % down
        for j=globminindex:-1:2
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(1:j));
            if(minindex == j-1 && ~isnan(minchi2) && (crumin-minchi2)>dchi2)         
                candidateindex_down = union(candidateindex_down, minindex); %R2013a compatible                
            end
            if(chi2s(j-1) < crumin-.01 && chi2s(j+1) < crumin-.01)
                if chi2s(j-1) < chi2s(j+1)
                    candidateindex_down = union(candidateindex_down, j-1);
                else
                    candidateindex_up = union(candidateindex_up, j+1);
                end
            end
        end
        
        d = diff(chi2s(1:globminindex));
        dsort = -sort(-d);  % biggest is 1st
        [~,inddmax] = max(d);
        if dsort(1) > dsort(2) + 5*(dsort(2)-dsort(3))  % if jump is 5x larger than last jump
            candidateindex_down = union(candidateindex_down,inddmax(1)-1);
        end

        % up
        for j=globminindex:1:(n-1)
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(j:n));
            if(minindex == 2 && ~isnan(minchi2) && (crumin-minchi2)>dchi2)
                candidateindex_up = union(candidateindex_up, (minindex+j-1)); %R2013a compatible
            end
        end
        d = diff(chi2s(n:-1:globminindex));
        dsort = -sort(-d);  % biggest is 1st
        [~,inddmax] = max(d);
        if dsort(1) > dsort(2) + 5*(dsort(2)-dsort(3))  % if jump is 5x larger than last jump
            candidateindex_up = union(candidateindex_up,n-inddmax(1)+1);
        end

        
        direction = [ones(size(candidateindex_down)) -ones(size(candidateindex_up))];
        candidateindex = [candidateindex_down candidateindex_up];
        trialP = [trialP; ar.ple.ps{jk}(candidateindex,:)]; 
        ar.ple.local_minima{jk} = trialP;

        if(~quick)
            if(size(trialP,1)>1)
                h2 = plePlot(jk, false, true, [], true, false);
                drawnow;
                trail_list = cell(1, size(trialP,1));
                for jp=1:size(trialP,1)
                    trail_list{jp} = sprintf('%s = %f', ar.ple.p_labels{jk}, trialP(jp,jk));
                end
                result = stringListChooser(trail_list, 1, true);
                candidateindex = candidateindex(result);
                trialP = trialP(result,:);
                direction = direction(result);
                close(h2);
            elseif(size(trialP,1)==0)
                fprintf('no trial point\n');
                return
            end
        else
            result = 1;
            candidateindex = candidateindex(result);
            trialP = trialP(result,:);
            direction = direction(result);
        end
    else
        % Find closest point
        [~,candidateindex]  = min( abs( ar.ple.ps{jk}(:,jk) - point ) );
        trialP              = ar.ple.ps{jk}(candidateindex,:);
        direction           = dr;
    end
    
    fprintf('trial point: %s = %f\n', ar.ple.p_labels{jk}, trialP(jk));
    
    if(direction<0)
        ntot = candidateindex;
    else
        ntot = n-candidateindex;
    end
    ncount = 0;
    arWaitbar(0);
    while(candidateindex+direction<=n && candidateindex+direction>0)
        ncount = ncount + 1;
        arWaitbar(ncount, ntot, sprintf('PLE#%i smoothing for %s', ...
                jk, strrep(ar.ple.p_labels{jk},'_', '\_')));
        pTrail = ar.ple.ps{jk}(candidateindex,:);
        pTrail(jk) = ar.ple.ps{jk}(candidateindex+direction,jk);
        
        feval(ar.ple.integrate_fkt, pTrail);
        [p, g] = feval(ar.ple.fit_fkt, jk);
        chi2 = feval(ar.ple.merit_fkt);
        
        if(chi2 < ar.ple.chi2s{jk}(candidateindex+direction))
            candidateindex = candidateindex+direction;
            ar.ple.ps{jk}(candidateindex,:) = p+0;
            ar.ple.chi2s{jk}(candidateindex) = chi2+0;
            ar.ple.gradient{jk}(candidateindex,:) = g+0;
            if(isfield(ar.ple,'violations'))
                ar.ple.chi2sviolations{jk}(candidateindex) = feval(ar.ple.violations);
            end
            if(isfield(ar.ple,'priors'))
                ar.ple.chi2spriors{jk}(candidateindex) = feval(ar.ple.priors, jk);
            end
            if(isfield(ar.ple,'priorsAll'))
                ar.ple.chi2spriorsAll{jk}(candidateindex) = feval(ar.ple.priorsAll);
            end
        else
            break;
        end
        
        if(ar.ple.showCalculation)
            plePlot(jk);
        end
    end
    arWaitbar(-1);
    
    fprintf('\n');
end

% reset parameters
feval(ar.ple.integrate_fkt, ar.ple.p);


% define new optimum
if(exist('chi2s','var') && ar.ple.merit-min(ar.ple.chi2s{jk}) > ar.ple.optimset_tol)
    [minchi2, iminchi2] = min(ar.ple.chi2s{jk});
    fprintf('PLE#%i found better optimum with chi^2 decrease of %e\n', jk, ...
        ar.ple.merit - minchi2);
    
    if(ar.ple.allowbetteroptimum)
        ar.ple.merit = minchi2;
        ar.ple.p = ar.ple.ps{jk}(iminchi2,:);
        feval(ar.ple.setoptim_fkt, ar.ple.ps{jk}(iminchi2,:));
    end
end

% save
pleSave(ar);

