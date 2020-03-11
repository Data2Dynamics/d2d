% pass = arCompareV(ar1,ar2,[silent])
% 
% compares ar.model(m).condition(c).vFineSimu in ar1 and ar2
%
%   silent - boolean for plotting (ar1 and sbmlexport in same plot) [false]
%   pass   - boolean, if relative difference <rtol, pass = 1 

function pass = arCompareV(ar1,ar2, silent)

% Tolerances for a pass
rtol = 1e-3;
atol = 1e-3;

if (nargin<3)
    silent = false;
end

maxDifference = [];
if size(ar1.model,2) ~= size(ar2.model,2)
    warning('arCompareV.m: Sizes of models not the same. Check model import.')
    pass = 0;
    return
end
for m=1:size(ar1.model,2)
    if size(ar1.model(m).condition,2) ~= size(ar2.model(m).condition,2)
        warning('arCompareV.m: Sizes of conditions not the same. Check model import.')
        pass = 0;
        return
    end
    for c=1:size(ar1.model(m).condition,2)
        t = ar2.model(m).condition(c).tFine;
        indcol = find(sum(ar2.model(m).condition(c).vFineSimu~=0,1)>0);

        if isempty(indcol)
            fprintf('arCompareV.m: No V found. Comparing V skipped.\n')
            pass = 0;
            return
        end

        subx = ceil(sqrt(length(indcol)-1));
        suby = ceil((length(indcol)-1)/subx);

        for i=2:length(ar2.model(m).v)
            if(length(ar2.model(m).v{i})==1)
                ar2.model(m).v{i} = [ar2.model(m).v{i},'_state'];
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
            sbmlsim = ar2.model(m).condition(c).vFineSimu(:,indcol(i));

            ind = strmatch(ar2.model(m).v{indcol(i)},ar1.model(m).v, 'exact'); %#ok
            if isempty(ind) 
                if (size(ar2.model(m).v,2) == size(ar1.model(m).v,2))
                    ind = indcol(i);
                    fprintf(['arCompareV.m: Names not consistent. Expecting ' ar2.model(m).v{indcol(i)} ' to be the same as ' ar1.model(m).v{ind} '. If not check your SBML export.\n']);
                else
                    if isempty(strmatch(ar2.model(m).v{indcol(i)},ar1.pLabel, 'exact')); %#ok
                        arFprintf(2, '%s from SBML export neither found as dynamic state nor as parameter.\n',ar2.model(m).v{indcol(i)})
                    end
                    continue
                end
            elseif(length(ind)>1)
                fprintf('%s from SBML export found multiple times.\n',ar2.model(m).v{indcol(i)})
            end
            d2dSim           = interp1(ar1.model(m).condition(c).tFine,ar1.model(m).condition(c).vFineSimu(:,ind),t);
            d2dFilt          = bsxfun(@max, d2dSim, atol);
            sbmlFilt    = bsxfun(@max, sbmlsim, atol);
            maxDifference(end+1) = max( ( (d2dFilt - sbmlFilt) ./ (sbmlFilt) ).^2 );

            if ( ~silent )
                subplot(subx,suby,i-1)
                set(gca,'FontSize',fs)
                plot(t, sbmlsim,'k');
                hold on
                plot(t, d2dSim,'r--')
                if dolegend ==1
                    set(legend(ar2.info.name,ar1.info.name),'FontSize',fs,'Interpreter','none');
                    dolegend = 0; % only once
                end
                title(strrep(ar2.model(m).v{indcol(i)},'_','\_'),'FontSize',fs);
                saveas(gcf,['arCompareVm' num2str(m) 'c' num2str(c)]);
            end
        end
    end
end

% Fit acceptable?
if ( max( maxDifference ) > rtol )
    pass = 0;
else
    pass = 1;
end




