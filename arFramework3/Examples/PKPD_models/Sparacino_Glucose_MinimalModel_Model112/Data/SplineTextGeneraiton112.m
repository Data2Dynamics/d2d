% open Simulated_magni_2006_data.csv, klick on "Import Selection"
% the spline is saved in SplineTextGeneration112.txt

n = length(Simulatedmagni2006data.TIME);
file = fopen('SplineTextGeneration112.txt','w');
fprintf(file,'number of knots: ' + string(n+2) + '\n\n');

fprintf(file,'TIME:\n[' + string(Simulatedmagni2006data.TIME(1)-1) + ', ');
for i = 1:n
    fprintf(file,string(Simulatedmagni2006data.TIME(i)) + ', ');
end
fprintf(file,string(Simulatedmagni2006data.TIME(n)+1) + ']\n\n');

fprintf(file,'INS:\n[' + string(Simulatedmagni2006data.INS(1)) + ', ');
for i = 1:n
    fprintf(file,string(Simulatedmagni2006data.INS(i)) + ', ');
end
fprintf(file,string(Simulatedmagni2006data.INS(n)) + ']\n\n');


