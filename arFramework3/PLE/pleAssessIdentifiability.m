% function for assessing whether parameters are point-wise identifiable
% while also checking that other parameters are not at the bounds
% 
% Used to set ar.ple.IDstatus_point_inBounds
% 
% INPUTS
% parsToCheck: cell array of parameter indices to check for identifiability [all parameters]
% alpha: confidence level for chi2 test [0.95]
% verbose: if true, print out which parameters are close to the bounds [false]
%

function [id,idR,idL] = pleAssessIdentifiability(parsToCheck,alpha,verbose)
global ar

if ~exist("alpha",'var') || isempty(alpha)
    alpha = 0.95;
end
if ~exist("parsToCheck",'var') || isempty(parsToCheck)
    parsToCheck =1:length(ar.ple.p);
end
if ~exist("verbose",'var') || isempty(verbose)
    verbose = false;
end

thresh = icdf('chi2',alpha,1);
closeVal = 0.001;

idL = NaN(size(ar.ple.chi2s));
idR = NaN(size(ar.ple.chi2s));
for i=1:length(ar.p)
    if ar.p(i) ~= 0
        [chi2min,imin] = min(ar.ple.chi2s{i});
        pi = ar.ple.ps{i}(imin(1),i);
    
        left = find(ar.ple.ps{i}(:,i)<pi);
        right = find(ar.ple.ps{i}(:,i)>pi);
        
        aboveL = ar.ple.chi2s{i}(left)'>(chi2min+thresh);
        aboveR = ar.ple.chi2s{i}(right)'>(chi2min+thresh);
        for j=parsToCheck{i}
            aboveL = aboveL & (ar.ple.ps{i}(left,j) < (ar.ple.ub(j)-closeVal));
            aboveL = aboveL & (ar.ple.ps{i}(left,j) > (ar.ple.lb(j)+closeVal));
            aboveR = aboveR & (ar.ple.ps{i}(right,j) < (ar.ple.ub(j)-closeVal));
            aboveR = aboveR & (ar.ple.ps{i}(right,j) > (ar.ple.lb(j)+closeVal));
            if ~verbose
                if(sum(ar.ple.ps{i}(right,j) > (ar.ple.ub(j)-closeVal)))
                    fprintf('Profile %i: parameter %i is close to upper bound (right)\n',i,j)
                end
                if(sum(ar.ple.ps{i}(right,j) < (ar.ple.lb(j)+closeVal)))
                    fprintf('Profile %i: parameter %i is close to lower bound (right)\n',i,j)
                end
                if(sum(ar.ple.ps{i}(left,j) > (ar.ple.ub(j)-closeVal)))
                    fprintf('Profile %i: parameter %i is close to upper bound (left)\n',i,j)
                end
                if(sum(ar.ple.ps{i}(left,j) < (ar.ple.lb(j)+closeVal)))
                    fprintf('Profile %i: parameter %i is close to lower bound (left)\n',i,j)
                end
            end
        end
        idL(i) = sum(aboveL)>0;
        idR(i) = sum(aboveR)>0;
    else
        idL(i) = 0;
        idR(i) = 0;
    end
end
id = idR & idL;
ar.ple.IDstatus_point_inBounds_alpha = alpha;
ar.ple.IDstatus_point_inBounds = id;



           