function vplSmooth
% vplSmooth
%
% Automatically smoothens peaks/jumps from the validation profile due to
% local otpima.
%
% Note an issue was detected where the prediction profile after smoothing
% may have a large space devoid of points. This is due to fitting the 
% validation profile from another local optimum, which changes the
% generated predictions. Proposal: Fill the prediction profile from the
% previous parameter region?

global ar

if(~isfield(ar,'vpl'))
    disp('Warning vplSmooth: Calculate validation profile first.')
    return
end

%Save the old ar-struct, to remove added data points later.
ar_old = arDeepCopy(ar);

try
    
    %Define shortcuts for refrence to general profile information:
    m = ar.vpl.general.m;
    d = ar.vpl.general.d;
    idpred = ar.vpl.general.idpred;
    tpred = ar.vpl.general.tpred;
    sigma = ar.vpl.general.sigma;
    
    %Set up some objects:
    dchi2 = 0.1;
    candidateindex_down = [];
    candidateindex_up = [];
    chi2s = ar.vpl.results.chi2;
    n = length(chi2s);
    [~, globminindex] = min(chi2s);
    
    % Find possible jump positions:
    
    % Downward direction:
    for j=globminindex:-1:2
        crumin = chi2s(j);
        [minchi2, minindex] = min(chi2s(1:j));
        % Possible Condition 0: Difference between current point and next
        %   point which is also the minimum on the set of points on the
        %   other side of the current point is too large
        % Remark: Why not just use the difference between current point
        %   and next point?
        if(minindex == j-1 && ~isnan(minchi2) && (crumin-minchi2)>dchi2)
            candidateindex_down = union(candidateindex_down, minindex); %R2013a compatible
        end
        % Possible Condition 1: Point is higher than left and right neighbour
        if(chi2s(j-1) < crumin-.005 && chi2s(j+1) < crumin-.005)
            if chi2s(j-1) < chi2s(j+1)
                candidateindex_down = union(candidateindex_down, j-1);
                %(downward) jump position is likely smaller than index j
            else
                candidateindex_up = union(candidateindex_up, j+1);
            end
        end
    end
    
    % Possible Condidition 2, properties:
    %   -Will in most cases find a jump in a normally unimodal profile
    %       if there is only one jump
    %   -Should in many cases find a jump in a normally multimodal profile,
    %       if there is only one jump
    %   -In the case of multiple jumps, even prominent jumps might not
    %       be found
    di = diff(chi2s(1:globminindex));
    dsort = -sort(-di);  % biggest is 1st
    [~,inddmax] = max(di);
    if dsort(1) > dsort(2) + 5*(dsort(2)-dsort(3))  % if jump is 5x larger than last jump
        candidateindex_down = union(candidateindex_down,inddmax(1)-1);
    end
    
    % upward direction:
    for j=globminindex:1:(n-1)
        crumin = chi2s(j);
        [minchi2, minindex] = min(chi2s(j:n));
        if(minindex == 2 && ~isnan(minchi2) && (crumin-minchi2)>dchi2)
            candidateindex_up = union(candidateindex_up, (minindex+j-1)); %R2013a compatible
        end
        if(chi2s(j-1) < crumin-.005 && chi2s(j+1) < crumin-.005)
            if chi2s(j-1) < chi2s(j+1)
                candidateindex_down = union(candidateindex_down, j-1);
            else
                candidateindex_up = union(candidateindex_up, j+1);
            end
        end
    end
    di = diff(chi2s(n:-1:globminindex));
    dsort = -sort(-di);  % biggest is 1st
    [~,inddmax] = max(di);
    if dsort(1) > dsort(2) + 5*(dsort(2)-dsort(3))  % if jump is 5x larger than last jump
        candidateindex_up = union(candidateindex_up,n-inddmax(1)+1);
    end
    
    %go up from candidateindex_down, go down from candidateindex_up
    direction = [ones(size(candidateindex_down)) -ones(size(candidateindex_up))];
    candidateindex = [candidateindex_down candidateindex_up];
    trialP = ar.vpl.results.ps(candidateindex,:);
    trialZ = ar.vpl.results.z(candidateindex);
    
    %Remove candidates which include a NaN parameter value
    q_notnan = (sum(isnan(trialP),2) == 0);
    direction_list = direction(q_notnan);
    candidateindex_list = candidateindex(q_notnan);
    
    %Add data point for modification later:
    arAddToData(m,d,idpred,tpred,0,2,1);
    iddata = length(ar.model(m).data(d).yExp(:,idpred));
    ar.model(m).data(d).yExpStd(iddata,idpred) = sigma;
    
    if isempty(candidateindex)
        disp('VPL smoothing SKIPPED: candidateindex empty.');
        ar = ar_old;
        return
    else
        %loop over all indices to automatically smoothen all jumps
        for ii = 1:length(candidateindex)
            candidateindex = candidateindex_list(ii);
            direction = direction_list(ii);
            
            fprintf('Trial point: z = %0.4g\n', trialZ(ii))
            
            while(candidateindex+direction<=n && candidateindex+direction>0)
                %Set up new data point and objective function:
                pTrail = ar.vpl.results.ps(candidateindex,:);
                z_new = ar.vpl.results.z(candidateindex + direction);
                ar.model(m).data(d).yExp(iddata,idpred) = z_new;
                                
                %Evaluate the step:
                ar.p = pTrail;
                arLink(true);
                arCalcMerit(ar.vpl.config.sensi,ar.p(ar.qFit==1)); %might be unnecessary
                arFit(true);
                arCalcMerit(ar.vpl.config.sensi,ar.p(ar.qFit==1));
                pred_new = ar.model(m).data(d).yExpSimu(iddata,idpred);
                chi2 = arGetMerit(true)-ar.vpl.general.chi2_norm;
                %New chi2 value normalized 
                
                if(chi2 < ar.vpl.results.chi2(candidateindex+direction))
                    fprintf('Found a better chi2 with decrease of %0.4g \n',...
                        ar.vpl.results.chi2(candidateindex+direction)-chi2);
                    candidateindex = candidateindex+direction;
                    ar.vpl.results.ps(candidateindex,:) = ar.p;
                    ar.vpl.results.chi2(candidateindex) = chi2;
                    ar.vpl.results.pred(candidateindex) = pred_new;
                    ar.vpl.results.ppl(candidateindex) = chi2-((z_new-pred_new)/sigma)^2;
                else
                    break;
                end
            end
        end
    end
catch exception
    fprintf('ERROR vplSmooth: Resetting ar struct. Error message: \n %s \n Line: %s \n' ,...
        exception.message, sprintf('%i, ',exception.stack.line));
    ar = ar_old;
    return
end
%Restore old struct with updated validation profile:
temp_struct = ar.vpl;
ar = ar_old;
ar.vpl = temp_struct;
%chi2dif is not yet updated:
n = ar.vpl.config.maxsteps;
chi2s =ar.vpl.results.chi2;
ar.vpl.test.chi2dif = ...
    [chi2s(1:n)-chi2s(2:(n+1));0;chi2s((n+2):end)-chi2s((n+1):(end-1))];
end



