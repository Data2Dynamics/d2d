% pass = arCompareX(ar1,ar2,[silent])
% 
% compares ar.model(m).condition(c).xFineSimu in ar1 and ar2
%
%   silent - boolean for plotting (ar1 and sbmlexport in same plot) [false]
%   pass   - boolean, if relative difference <rtol, pass = 1 

function pass = arCompareX(ar1,ar2, silent)

% Tolerances for a pass
rtol = 1e-3;
atol = 1e-4;

if (nargin<3)
    silent = false;
end

maxDifference = [];
if size(ar1.model,2) ~= size(ar2.model,2)
    warning('arCompareX.m: Sizes of models not the same. Check model import.')
    pass = 0;
    return
end
for m=1:size(ar1.model,2)
    if size(ar1.model(m).condition,2) ~= size(ar2.model(m).condition,2)
        warning('arCompareX.m: Sizes of conditions not the same. Check model import.')
        pass = 0;
        return
    end
    for c=1:size(ar1.model(m).condition,2)
        t = ar2.model(m).condition(c).tFine;
        indcol = find(sum(ar2.model(m).condition(c).xFineSimu~=0,1)>0);

        if isempty(indcol)
            warning('arCompareX.m: No X found. Comparing X skipped.')
            pass = 0;
            return
        end

        subx = ceil(sqrt(length(indcol)-1));
        suby = ceil((length(indcol)-1)/subx);

        for i=2:length(ar2.model(m).x)
            if(length(ar2.model(m).x{i})==1)
                ar2.model(m).x{i} = [ar2.model(m).x{i},'_state'];
            end
        end

        if(length(indcol)<16)
            fs = 10;
        elseif(length(indcol)<25)
            fs = 8;
        elseif(length(indcol)<36)
            fs = 7;
        else
            fs = 6;
        end


        %%
        close all
        dolegend = 1;
        for i=length(indcol):-1:2
            sbmlsim = ar2.model(m).condition(c).xFineSimu(:,indcol(i));
            if ( ~silent )
                subplot(subx,suby,i-1)
                set(gca,'FontSize',fs)
                plot(t, sbmlsim,'k');
                hold on
            end

            ind = strmatch(ar2.model(m).x{indcol(i)},ar1.model(m).x, 'exact'); %#ok
            if isempty(ind) 
                if (size(ar2.model(m).x,2) == size(ar1.model(m).x,2))
                    ind = indcol(i);
                    warning(['arCompareX.m: Names not consistent. Expecting ' ar2.model(m).x{indcol(i)} ' to be the same as ' ar1.model(m).x{ind} '. If not check your SBML export.']);
                else
                    if isempty(strmatch(ar2.model(m).x{indcol(i)},ar1.pLabel, 'exact')); %#ok
                        arFprintf(2, '%s from SBML export neither found as dynamic state nor as parameter.\n',ar2.model(m).x{indcol(i)})
                    end
                    continue
                end
            elseif(length(ind)>1)
                warning('%s from SBML export found multiple times.\n',ar2.model(m).x{indcol(i)})
            end
            d2dSim           = interp1(ar1.model(m).condition(c).tFine,ar1.model(m).condition(c).xFineSimu(:,ind),t);
            d2dFilt          = bsxfun(@max, d2dSim, atol);
            sbmlFilt    = bsxfun(@max, sbmlsim, atol);
            maxDifference(end+1) = max( ( (d2dFilt - sbmlFilt) ./ (sbmlFilt) ).^2 );

            if ( ~silent )
                plot(t, d2dSim,'r--')
                if dolegend ==1
                    set(legend(ar2.info.name,ar1.info.name),'FontSize',fs,'Interpreter','none');
                    dolegend = 0; % only once
                end
            end

            if ( ~silent )
                xlim([0,100]);
                title(strrep(ar2.model(m).x{indcol(i)},'_','\_'),'FontSize',fs);
            end
        end
        if ( ~silent )
            saveas(gcf,['arCompareXm' num2str(m) 'c' num2str(c)]);
        end
    end
end

% Fit acceptable?
if ( max( maxDifference ) > rtol )
    pass = 0;
else
    pass = 1;
end




