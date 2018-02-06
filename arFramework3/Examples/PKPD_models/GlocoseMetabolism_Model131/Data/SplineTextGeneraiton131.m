% open SimuInsulinAction.csv, klick on "Import Selection"
% the spline is saved in SplineTextGeneration131.txt

n = length(SimuInsulinAction.TIME);
file = fopen('SplineTextGeneration131.txt','w');
fprintf(file,'number of knots: ' + string(n+2) + '\n\n');

fprintf(file,'TIME:\n[' + string(SimuInsulinAction.TIME(1)-1) + ', ');
for i = 1:n
    fprintf(file,string(SimuInsulinAction.TIME(i)) + ', ');
end
fprintf(file,string(SimuInsulinAction.TIME(n)+1) + ']\n\n');

fprintf(file,'INS:\n[' + string(SimuInsulinAction.INS(1)) + ', ');
for i = 1:n
    fprintf(file,string(SimuInsulinAction.INS(i)) + ', ');
end
fprintf(file,string(SimuInsulinAction.INS(n)) + ']\n\n');

fprintf(file,'Rinf:\n[' + string(SimuInsulinAction.Rinf(1)) + ', ');
for i = 1:n
    fprintf(file,string(SimuInsulinAction.Rinf(i)) + ', ');
end
fprintf(file,string(SimuInsulinAction.Rinf(n)) + ']\n\n');
