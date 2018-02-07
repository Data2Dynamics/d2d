% open data127original.csv and click "Import Selection"

n = length(data127original.TIME);
listID = zeros(n,2);

% create datapoints for the individual parameters
for i = 1:length(listID)
    listID(i,1) = data127original.ID(i);
    listID(i,2) = NaN;
end

id = listID(1,1)-1;
for j = (1:n)
    if listID(j,1) > id
        listID(j,2) = 0;
        id = listID(j,1);
    end
end

% put everything together in one table   
data127original.p_ID_Kgrw1 = listID(:,1);
data127original.p_ID_TS0 = listID(:,1);
data127original.p_ID_AMT = listID(:,1);
data127original.mean_p_ID_Kgrw1 = listID(:,2);
data127original.mean_p_ID_TS0 = listID(:,2);
data127original.mean_p_ID_AMT = listID(:,2);

% change "0" to "NaN" in DV for MDV = 1
for i = 1:n
    if data127original.MDV(i) == 1
        data127original.DV(i) = NaN;
    end
end

% export
writetable(data127original,'data127.csv')

% print a list for AMT (individual initial condition for dep)
list = data127original.AMT;
list = list(list>0);
join(string(list),', ')