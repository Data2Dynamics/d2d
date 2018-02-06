% open Simulated_data_PK_reduced.csv and click "Import Selection"

n = length(SimulateddataPKreduced.TIME);
tableExp = SimulateddataPKreduced;
listId = zeros(n,4);

for i = 1:n
    listId(i,1) = tableExp.ID(i);
    listId(i,2) = NaN;
    listId(i,3) = NaN;
    listId(i,4) = NaN;
end

id = listId(1,1)-1;
for i = (1:n)
    if listId(i,1) > id
        listId(i,2) = 0;
        listId(i,3) = tableExp.DDI(i);
        listId(i,4) = tableExp.CGCL(i);
        id = listId(i,1);
    end
end

% change order of columns that the first one is time (or TALD, time after
% liquid dose)
tableExp = [tableExp(:,9), tableExp(:,1:(end-1))];

% add the columns for the independant parameters to the table
tableExp.p_ID_CL = listId(:,1);
tableExp.p_ID_V2 = listId(:,1);
tableExp.p_ID_Q = listId(:,1);
tableExp.p_ID_DDIxxx = listId(:,1);
tableExp.p_ID_CGCLxxx = listId(:,1);

tableExp.mean_p_ID_CL = listId(:,2);
tableExp.mean_p_ID_V2 = listId(:,2);
tableExp.mean_p_ID_Q = listId(:,2);
tableExp.mean_p_ID_DDIxxx = listId(:,2);
tableExp.mean_p_ID_CGCLxxx = listId(:,2);


% get a list for the DDI and CGCL parameters
list = listId(:,3);
list = list(list>=0);
join(string(list),', ')
list = listId(:,4);
list = list(list>=0);
join(string(list),', ')


% check if DDI and CGCL are constant for each Id (they are)
m1 = tableExp.ID(1);
m2 = tableExp.ID(end);
for i = m1:m2
    ind = find(tableExp.ID == i);
    if max(tableExp.DDI(ind)) - min(tableExp.DDI(ind)) > 10^-9
        disp('ERROR DDI: ' + string(i));
    end
    if max(tableExp.CGCL(ind)) - min(tableExp.CGCL(ind)) > 10^-9
        disp('ERROR CGCL: ' + string(i));
    end
end

writetable(tableExp,'data118.csv')