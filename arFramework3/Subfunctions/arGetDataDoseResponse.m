function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, data_qFit, yExpHl] = ...
    arGetDataDoseResponse(jm, ds, ttime, dLink, logplotting_xaxis, jtype)
global ar

zero_break = [];
data_qFit = true;

if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

ccount = 1;
for jd = ds
    
    % get condition jcondi
    data_qFit = data_qFit & ar.model(jm).data(jd).qFit;
    qt = ar.model(jm).data(jd).tExp == ttime;
    for jc = 1:length(ar.model(jm).data(jd).condition)
        if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ...
                ar.model(jm).data(jd).response_parameter))
            jcondi = jc;
        end
    end
    
    jc = ar.model(jm).data(jd).cLink;
    
    for jt = find(qt')
        
        % collect x axis values
        if(logplotting_xaxis)
            t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(jcondi).value)); %#ok<AGROW>
        else
            t(ccount,1) = str2double(ar.model(jm).data(jd).condition(jcondi).value); %#ok<AGROW>
        end
        
        % calculatte zero break
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = dLink
                if(logplotting_xaxis)
                    if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value))))
                        doses(end+1) = ...
                            log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value)); %#ok<AGROW>
                    end
                else
                    doses(end+1) = str2double(ar.model(jm).data(jd2).condition(jcondi).value); %#ok<AGROW>
                end
            end
            doses = unique(doses); %R2013a compatible
            if(length(doses)>1)
                t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            else
                t(ccount,1) = doses(1) - 0.1; %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            end
        end
        tExp(ccount,1) = t(ccount,1); %#ok<AGROW>
        
        % trajectories and error bands
        if(jtype==1)
            [tdiffmin, jtfine] = min(abs(ar.model(jm).data(jd).tFine - ...
                ar.model(jm).data(jd).tExp(jt))); %#ok<ASGLU>
            y(ccount,:) = ar.model(jm).data(jd).yFineSimu(jtfine,:); %#ok<AGROW>
            ystd(ccount,:) = ar.model(jm).data(jd).ystdFineSimu(jtfine,:); %#ok<AGROW>
        elseif(jtype==2)
            [tdiffmin, jtfine] = min(abs(ar.model(jm).condition(jc).tFine - ...
                ar.model(jm).data(jd).tExp(jt))); %#ok<ASGLU>
            y(ccount,:) = [ar.model(jm).condition(jc).uFineSimu(jtfine,:) ...
                ar.model(jm).condition(jc).xFineSimu(jtfine,:) ...
                ar.model(jm).condition(jc).zFineSimu(jtfine,:)]; %#ok<AGROW>
            ystd = [];
            
        elseif(jtype==3)
            [tdiffmin, jtfine] = min(abs(ar.model(jm).condition(jc).tFine - ...
                ar.model(jm).data(jd).tExp(jt))); %#ok<ASGLU>
            y(ccount,:) = [ar.model(jm).condition(jc).vFineSimu(jtfine,:)]; %#ok<AGROW>
            ystd = [];
        end

        % data
        if(jtype==1)
            yExp(ccount,:) = ar.model(jm).data(jd).yExp(jt,:); %#ok<AGROW>
            yExpHl(ccount,:) = NaN; %#ok<AGROW>
            if(isfield(ar.model(jm).data(jd),'highlight'))
                if(ar.model(jm).data(jd).highlight(jt,:)~=0)
                    yExpHl(ccount,:) = yExp(ccount,:);                 %#ok<AGROW>
                end
            end
            if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==-1) || ...
                        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==-1) )
                yExpStd(ccount,:) = ar.model(jm).data(jd).yExpStd(jt,:); %#ok<AGROW>
            else
                yExpStd(ccount,:) = ar.model(jm).data(jd).ystdExpSimu(jt,:); %#ok<AGROW>
            end
        else
            tExp = [];
            yExp = [];
            yExpStd = [];
            yExpHl = [];
        end
        
        % confidence bands
        if(jtype==1)
            if(isfield(ar.model(jm).data(jd), 'yExpUB'))
                lb(ccount,:) = ar.model(jm).data(jd).yFineLB(jtfine,:); %#ok<AGROW>
                ub(ccount,:) = ar.model(jm).data(jd).yFineUB(jtfine,:); %#ok<AGROW>
            else
                lb = [];
                ub = [];
            end
            
        elseif(jtype==2)
            if(isfield(ar.model(jm).data(jd), 'xExpUB'))
                lb(ccount,:) = [ar.model(jm).condition(jc).uFineLB(jtfine,:) ...
                    ar.model(jm).condition(jc).xFineLB(jtfine,:) ...
                    ar.model(jm).condition(jc).zFineLB(jtfine,:)]; %#ok<AGROW>
                ub(ccount,:) = [ar.model(jm).condition(jc).uFineUB(jtfine,:) ...
                    ar.model(jm).condition(jc).xFineUB(jtfine,:) ...
                    ar.model(jm).condition(jc).zFineUB(jtfine,:)]; %#ok<AGROW>
            else
                lb = [];
                ub = [];
            end
            
        elseif(jtype==3)
            if(isfield(ar.model(jm).data(jd), 'vExpUB'))
                lb(ccount,:) = ar.model(jm).condition(jc).vFineLB(jtfine,:); %#ok<AGROW>
                ub(ccount,:) = ar.model(jm).condition(jc).vFineUB(jtfine,:); %#ok<AGROW>
            else
                lb = [];
                ub = [];
            end
            
        end
        
        ccount = ccount + 1;
    end
end

if(~isempty(yExp))
    [tExp,itexp] = sort(tExp);
    yExp = yExp(itexp,:);
    yExpHl = yExpHl(itexp,:);
else
    data_qFit = [];
end

[t,it] = sort(t);
y = y(it,:);
if(~isempty(ystd))
    ystd = ystd(it,:);
end
if(~isempty(lb))
    lb = lb(it,:);
    ub = ub(it,:);
end