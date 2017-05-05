% Smoothen jumps to local optima
%
% Usage: 
%   function arPLESmooth(jk, quick, point, dr)
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

function varargout = arPLESmooth(varargin)

global ar

if(nargin==0)
    qglobalar = true;
    jk = find(ar.qFit==1);
    quick = false;
else
    if(isstruct(varargin{1}))
        qglobalar = false;
        ar = varargin{1};
        if(nargin==1)
            jk = find(ar.qFit==1);
            quick = false;
        elseif(nargin==2)
            jk = varargin{2};
            quick = false;
        elseif(nargin==3)
            jk = varargin{2};
            quick = varargin{3};
        end
    else
        qglobalar = true;
        if(nargin==1)
            jk = varargin{1};
            quick = false;
        elseif(nargin==2)
            jk = varargin{1};
            quick = varargin{2};
        end
    end
end

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 
if(~exist('quick','var'))
    quick = false;
end
if(nargin<1)
    fprintf('PLE smoothing for %i parameters ...\n', sum(ar.qFit))
    jindex = find(ar.qFit);
    for j=1:length(jindex)
        if(qglobalar)
            arPLESmooth(jindex(j));
        else
            ar = arPLESmooth(jindex(j));
        end
    end
    return
elseif(length(jk)>1)
    fprintf('PLE smoothing for %i parameters ...\n', length(jk))
    for j=1:length(jk)
        if(qglobalar)
            arPLESmooth(jk(j));
        else
            ar = arPLESmooth(jk(j));
        end
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
    fprintf('PLE#%i smoothing SKIPPED: parameter %s is fixed\n', jk, ar.pLabel{jk});
    return;
elseif(jk <= length(ar.ple.chi2s) && ~isempty(ar.ple.chi2s{jk}))
    fprintf('PLE#%i smoothing for parameter %s\n', jk, ar.pLabel{jk});
        
    trialP = zeros(0,length(ar.p));

    n = length(ar.ple.chi2s{jk});
    chi2s = ar.ple.chi2s{jk};
    [~, globminindex] = min(chi2s);
    
    hasTrialpoint = true;
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
        trialP = [trialP; ar.ple.ps{jk}(candidateindex,:)]; 
        ar.ple.local_minima{jk} = trialP;

        if(~quick)
            if(size(trialP,1)>1)
                trail_list = cell(1, size(trialP,1));
                for jp=1:size(trialP,1)
                    trail_list{jp} = sprintf('PLE#%i %s = %f', jk, ar.pLabel{jk}, trialP(jp,jk));
                end
                result = stringListChooser(trail_list, 1, true);
                candidateindex = candidateindex(result);
                trialP = trialP(result,:);
                direction = direction(result);
            elseif(size(trialP,1)==0)
                fprintf('PLE#%i no trial point\n', jk);
                hasTrialpoint = false;
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
    
    if(hasTrialpoint)
        fprintf('PLE#%i trial point: %s = %f\n', jk, ar.pLabel{jk}, trialP(jk));
        
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
                jk, strrep(ar.pLabel{jk},'_', '\_')));
            pTrail = ar.ple.ps{jk}(candidateindex,:);
            pTrail(jk) = ar.ple.ps{jk}(candidateindex+direction,jk);
            
            ar = arCalcMerit(ar, true, pTrail(ar.qFit==1));
            ar = arFit(ar, true);
            p = ar.p;
            chi2 = arGetMerit;
            
            if(chi2 < ar.ple.chi2s{jk}(candidateindex+direction))
                candidateindex = candidateindex+direction;
                ar.ple.ps{jk}(candidateindex,:) = p+0;
                ar.ple.chi2s{jk}(candidateindex) = chi2+0;
            else
                break;
            end
        end
        arWaitbar(-1);
    end
end

% reset parameters
ar = arCalcMerit(ar, true, ar.ple.pStart(ar.qFit==1));

% define new optimum
if(exist('chi2s','var') && arGetMerit(true)-min(ar.ple.chi2s{jk}) > 1e-1)
    minchi2 = min(ar.ple.chi2s{jk});
    fprintf('PLE#%i found better optimum with chi^2 decrease of %e\n', jk, ...
        ar.ple.chi2Reset - minchi2);
end

if(nargout>0 && ~qglobalar)
    varargout{1} = ar;
else
    varargout = {};
end


