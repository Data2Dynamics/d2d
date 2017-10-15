
y_tmp = ar.model(1).data(1).yExp(4,1);
ys = linspace(y_tmp-2,y_tmp+2,101);
for i=1:101
    ar.model(1).data(1).yExp(4,1) = ys(i);
    arCalcMerit
    chi2s(i) = ar.chi2;
end
ar.model(1).data(1).yExp(4,1) = y_tmp;
plot(ys,chi2s)
set(gcf,'Color','w')
xlabel('data point')
ylabel('chi2 value')
title('check of Laplace distribution')