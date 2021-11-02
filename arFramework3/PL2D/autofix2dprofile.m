function autofix2dprofile(fixmodes)

% autofix2dprofile(fixmodes)
%
% Auto-fixes all deficiencies of 2d profile.
%
% fixmodes: Determines which of the three main fixing routines is employed.
%               [Default: logical([1,1,1])]
%           Mode 1: rmProfileJumps2d
%               Removes discontinuities occuring between whole profiles by
%               recalculating these profiles. Discontinuties are found
%               automatically by supplementary functions 
%               SmoothingPoints2d_Excess, SmoothingPoints2d_Bounds and
%               SmoothingPoints2d_Inside. This is the main correction 
%               function.
%           Mode 2: extend2dprofile
%               Extends profiles which did not converge to either the
%               confidence threshold or parameter bounds.
%           Mode 3: rmInsideJumps2d
%               Checks if there are any discontinuites due to non-convergred 
%               fits for interior 2d profile points and corrects them.
%
% It is recommended to run this function for every calculated 2d-profile.
% If no correction is needed, this function will skip all corrections.
% However, searching for discontinuities still takes some time even if there
% are none.
%
% Mode 1 is complicated to use as a separate function call, but extending
% profiles with extend2dprofiles or removing jumps due to non-converged fits
% with rmInsideJumps2d is easy enough.
%
% See also: gen2d, rmProfileJumps2d, extend2dprofile, rmInsideJumps2d

global ar

if (~isfield(ar,'ple2d')) || (~isfield(ar.ple2d,'raw'))
    fprintf('\n ERROR autofix2dprofile: Generate a 2d profile first. \n')
end
if ~exist('fixmodes','var') || isempty(fixmodes)
    fixmodes = logical([1,1,1]);
end

xnodes_reset = ar.ple2d.config.plot.xnodes;
ar.ple2d.config.plot.xnodes = 2000;

if fixmodes(1)
    
    fprintf('\n Removing discontinuities between profiles ... \n')
    
    for lbub = [-1,1]
        if lbub == 1
            txt = 'upper';
        elseif lbub == -1
            txt = 'lower';
        end
        
        iter = 0;
        while iter < ar.ple2d.config.autofix.niter + 1
            
            q1 = 0;
            q2 = 0;
            
            if iter > ar.ple2d.config.autofix.niter - 1
                fprintf([' \n WARNING: autofix2dprofile: Outside Profile',...
                    ' smoothing for %s bound did not finish in %i',...
                    ' iterations. \n'],txt,ar.ple2d.config.autofix.niter)
                break
            end
            iter = iter +1;
            
            % Try to remove profile jumps with excess method:
            [~,indys,dirs] = SmoothingPoints2d_Excess(lbub,1);
            for ii = 1:length(indys)
                fprintf(['\n Attempt %s bound profile smoothing for',...
                    ' trial point number %i in iteration %i',...
                    ' (excess check)... \n'],txt,ii,iter)
                rmProfileJumps2d(1,indys(ii),dirs(ii),lbub);
            end
            % Check whether same profile jumps are found. If yes,
            % no further correction was performed.
            [~,indys_new] = SmoothingPoints2d_Excess(lbub,0);
            if ((length(indys) == length(indys_new)) && ...
                    (sum(indys == indys_new) == length(indys))) || ...
                    isempty(indys_new)
                q1 = 1;
            end
            
            % Try to remove profile jumps with inside method:
            [~,indys,dirs] = SmoothingPoints2d_Inside(1);
            for ii = 1:length(indys)
                fprintf(['\n Attempt %s bound profile smoothing for',...
                    ' trial point number %i in iteration %i',...
                    ' (inside check)... \n'],txt,ii,iter)
                rmProfileJumps2d(3,indys(ii),dirs(ii),lbub);
            end
            [~,indys_new] = SmoothingPoints2d_Inside(0);
            if ((length(indys) == length(indys_new)) && ...
                    (sum(indys == indys_new) == length(indys))) || ...
                    isempty(indys_new)
                q2 = 1;
            end
            
            % If neither excess nor inside method improves the 2dprofile
            % further, switch to bound method:
            if (q1 == 1) && (q2 == 1)
                
                % Finds bound method smoothing points (computationally
                % comparatively expensive)
                [~,indys,dirs] = SmoothingPoints2d_Bounds(lbub,1);
                
                % This check is enough since the method proposes only
                % actually confirmed discontinuities.
                if isempty(indys)
                    fprintf('\n No more %s bound profile jumps found \n',txt)
                    break
                end
                
                % Removes profile jumps with bound method:
                for ii = 1:length(indys)
                    fprintf(['\n Attempt %s bound profile smoothing for',...
                        ' trial point number %i in iteration %i',...
                        ' (bounds check)... \n'],txt,ii,iter)
                    rmProfileJumps2d(2,indys(ii),dirs(ii),lbub);
                end
                
                % Now the next iteration starts. Note that the structure is
                % hierarchical, excess and inside method is preferred over
                % bound method, because finding the discontinuities is much
                % faster in the former methods, although the discotinutities
                % they find need not completely coincide.
            end
            
        end
        
    end
    
end

if fixmodes(2)
    
    fprintf('\n Extend non-converged profiles ... \n')
    
    extend2dprofile;
    
end

if fixmodes(3)
    
    fprintf('\n Removing inside discontinuities ... \n')
    
    % This function iterates until all discontinuities to the set threshold
    % are fixed. Efficiency can be considerably improved if the
    % discontinuites which were not fixed in the last iteration are not
    % checked again in the next iteration.
    
    iter = 0;
    count_improvement = 0;
    count_improvement_tmp = 1;
    while (count_improvement_tmp ~= 0)
        
        if iter > ar.ple2d.config.autofix.niter - 1
            fprintf([' \n WARNING: autofix2dprofile: Removing single',...
                ' discontinuities did not finish in %i',...
                ' iterations. \n'],txt,ar.ple2d.config.autofix.niter)
            break
        end
        iter = iter +1;
        
        [count_improvement_tmp,~] = rmInsideJumps2d;
        
        fprintf('\n %i discontinuities were corrected in iteration %i. \n',...
            count_improvement_tmp,iter)
        
        count_improvement = count_improvement_tmp + count_improvement;
    end
    fprintf('\n %i discontinuities were corrected overall. \n',count_improvement);
    
end

if fixmodes(2)
    
    fprintf('\n Extend non-converged profiles ... \n')
    
    extend2dprofile;
    
end

ar.ple2d.config.plot.xnodes = xnodes_reset;
smooth2d;

end
