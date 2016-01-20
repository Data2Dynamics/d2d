% [units,erg] = arCalcParameterUnitsViaSymmetries(m)
%
% Example:
% [units,erg] = arCalcParameterUnitsViaSymmetries;
% feval(erg.showResult)

function [units,erg] = arCalcParameterUnitsViaSymmetries(m)
if(~exist('m','var') || isempty(m))
    m = 1;
end

global ar

filenames = arExportModelToDaniel({'model','obsCalcUnits'});

warning(' The tool currently is only able to calculate the concentration units. Time units are negleted. Parameters without concentration dimension do not occur in the output');


python_string = sprintf('python symmetryDetection.py %s %s',filenames{1},filenames{2});
disp(python_string)
system(python_string);

fileResult = [filenames{1},'_result.txt'];
[group,variable,trsf] = arReadSymmetryDetectionOutput(fileResult);

%% edit output von symmetry tool
trsf2 = strrep(trsf,'*exp(','^(');
trsf2 = strrep(trsf2,'.0*','*');

units = cell(size(trsf2));
for i=1:length(trsf2)
    units{i} = subs(sym(trsf2{i}),'epsilon',1);
end

%% determine conc. units within each group, check if conc_units are inverse.
guni = unique(group);  % group names

conc_units = regexp(variable,'conc_unit\d+','match');  % unit labels
conc_units = [conc_units{find(~cellfun(@isempty,conc_units))}];  % remove empty

doInverse = zeros(size(units));

unitGroup = NaN(size(units));  % Indicates whether and which group of concentration units
for g=1:length(guni)
    inGroup = find(ismember(group,guni{g}));  % indices belonging to the group g
    mom_conc_unit = intersect(conc_units,variable(inGroup)); % current unit labels
    
    if(~isempty(mom_conc_unit))  % if the group corresponds by fixing a concentration unit
        ind =  strmatch(mom_conc_unit, variable(inGroup),'exact');    % which variable is the concentration unit label
        if(units{inGroup(ind)} == sym([variable{inGroup(ind)},'^1']))
            doInverse(inGroup) = 1;
        end
        
        for i=1:length(inGroup)
            unitGroup(inGroup(i)) = g;
        end
        
        [~,ia] = intersect(ar.model(m).x,variable(inGroup));
        
        currentUnit = unique(ar.model(m).xUnits(ia,1));
        if iscell(currentUnit) && length(currentUnit)>1
            error('Units in  ar.model.xUnits(:,1) are not unique for %s ',sprintf('%s ',variable(inGroup)));
        end
        
        for i=1:length(inGroup)
            units{inGroup(i)} = subs(units{inGroup(i)},variable{inGroup(i)},currentUnit);
        end
        
    else % syms without relationship to conc. units
        for i=1:length(inGroup)
            units{inGroup(i)} = sym(1);
        end
    end
end

doInv = find(doInverse);
for i=1:length(doInv)
    units{doInv(i)} = 1/units{doInv(i)};
end


%% print result:
    function showResult
        guni = unique(erg.group);
        for gg=1:length(guni)
            ind = strmatch(guni{gg},erg.group,'exact');
            fprintf('\nGroup %s:\n',guni{gg});
            for ii=1:length(ind)
                if(isempty(intersect(erg.variable{ind(ii)},erg.conc_units)))
                    fprintf('[%s]\t=\t%s\n',erg.variable{ind(ii)},char(erg.units{ind(ii)}));
                end
            end
        end
    end


%% make result struct
erg.variable = variable;
erg.trsf = trsf;
erg.trsf2 = trsf2;
erg.units = units;
erg.unit_str = cellfun(@char,erg.units,'UniformOutput',false);
erg.group = group;
erg.conc_units = conc_units;
erg.unitGroup = unitGroup;
erg.showResult = @showResult;
feval(erg.showResult)

end