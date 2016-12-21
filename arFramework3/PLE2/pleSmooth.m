% Smoothen jumps to local optima
%
% Usage: 
%   function pleSmooth(jk, quick, point, dr)
% 
% Mandatory parameter:
%   jk         parameter to smoothen
%
% Optional parameters
%   quick      skip selection and pick the first trial point (default = 0) 
%   point      provide custom parameter value to start at (closest point
%              will be chosen)
%   dr         initial search direction (-1 or 1; required when point is 
%              specified)

function pleSmooth(jk, quick, point, dr)

global pleGlobals;

if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 
if(~isfield(pleGlobals, 'showCalculation'))
    pleGlobals.showCalculation = true;
end
if(~exist('quick','var'))
    quick = false;
end
if(nargin<1)
    fprintf('PLE smoothing for %i parameters ...\n', sum(pleGlobals.q_fit))
    jindex = find(pleGlobals.q_fit);
    for j=1:length(pleGlobals.chi2s)
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

if(~pleGlobals.q_fit(jk))
    fprintf('\nPLE#%i smoothing SKIPPED: parameter %s is fixed\n', jk, pleGlobals.p_labels{jk});
    return;
elseif(jk <= length(pleGlobals.chi2s) && ~isempty(pleGlobals.chi2s{jk}) && pleGlobals.IDstatus(jk)~=4)
    fprintf('\nPLE#%i smoothing for parameter %s\n', jk, pleGlobals.p_labels{jk});
    
    
    trialP = zeros(0,length(pleGlobals.p));

    n = length(pleGlobals.chi2s{jk});
    chi2s = pleGlobals.chi2s{jk};
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

        % up
        for j=globminindex:1:(n-1)
            crumin = chi2s(j);
            [minchi2, minindex] = min(chi2s(j:n));
            if(minindex == 2 && ~isnan(minchi2) && (crumin-minchi2)>dchi2)
                candidateindex_up = union(candidateindex_up, (minindex+j-1)); %R2013a compatible
            end
        end

        direction = [ones(size(candidateindex_down)) -ones(size(candidateindex_up))];
        candidateindex = [candidateindex_down candidateindex_up];
        trialP = [trialP; pleGlobals.ps{jk}(candidateindex,:)]; 
        pleGlobals.local_minima{jk} = trialP;

        if(~quick)
            if(size(trialP,1)>1)
                h2 = plePlot(jk, false, true, [], true, false);
                drawnow;
                trail_list = cell(1, size(trialP,1));
                for jp=1:size(trialP,1)
                    trail_list{jp} = sprintf('%s = %f', pleGlobals.p_labels{jk}, trialP(jp,jk));
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
        [~,candidateindex]  = min( abs( pleGlobals.ps{jk}(:,jk) - point ) );
        trialP              = pleGlobals.ps{jk}(candidateindex,:);
        direction           = dr;
    end
    
    fprintf('trial point: %s = %f\n', pleGlobals.p_labels{jk}, trialP(jk));
    
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
                jk, strrep(pleGlobals.p_labels{jk},'_', '\_')));
        pTrail = pleGlobals.ps{jk}(candidateindex,:);
        pTrail(jk) = pleGlobals.ps{jk}(candidateindex+direction,jk);
        
        feval(pleGlobals.integrate_fkt, pTrail);
        [p, g] = feval(pleGlobals.fit_fkt, jk);
        chi2 = feval(pleGlobals.merit_fkt);
        
        if(chi2 < pleGlobals.chi2s{jk}(candidateindex+direction))
            candidateindex = candidateindex+direction;
            pleGlobals.ps{jk}(candidateindex,:) = p+0;
            pleGlobals.chi2s{jk}(candidateindex) = chi2+0;
            pleGlobals.gradient{jk}(candidateindex,:) = g+0;
            if(isfield(pleGlobals,'violations'))
                pleGlobals.chi2sviolations{jk}(candidateindex) = feval(pleGlobals.violations);
            end
            if(isfield(pleGlobals,'priors'))
                pleGlobals.chi2spriors{jk}(candidateindex) = feval(pleGlobals.priors, jk);
            end
            if(isfield(pleGlobals,'priorsAll'))
                pleGlobals.chi2spriorsAll{jk}(candidateindex) = feval(pleGlobals.priorsAll);
            end
        else
            break;
        end
        
        if(pleGlobals.showCalculation)
            plePlot(jk);
        end
    end
    arWaitbar(-1);
    
    fprintf('\n');
end

% reset parameters
feval(pleGlobals.integrate_fkt, pleGlobals.p);


% define new optimum
if(exist('chi2s','var') && pleGlobals.chi2-min(pleGlobals.chi2s{jk}) > pleGlobals.optimset_tol)
    [minchi2, iminchi2] = min(pleGlobals.chi2s{jk});
    fprintf('PLE#%i found better optimum with chi^2 decrease of %e\n', jk, ...
        pleGlobals.chi2 - minchi2);
    
    if(pleGlobals.allowbetteroptimum)
        pleGlobals.chi2 = minchi2;
        pleGlobals.p = pleGlobals.ps{jk}(iminchi2,:);
        feval(pleGlobals.setoptim_fkt, pleGlobals.ps{jk}(iminchi2,:));
    end
end

% save
if(~exist([cd '/' pleGlobals.savePath], 'dir'))
    mkdir([cd '/' pleGlobals.savePath])
end
doSave(pleGlobals);

%% don't use the global pleGlobals in this function to prevent that manipulations have effect on the global variable
function    doSave(pleGlobals)
pleGlobals.fighandel_multi = [];    % remove handle to
save([pleGlobals.savePath '/results.mat'], 'pleGlobals');
