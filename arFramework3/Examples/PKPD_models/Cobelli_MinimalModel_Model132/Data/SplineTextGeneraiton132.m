% open Real_minimalModel.csv, klick on "Import Selection"
% the spline is saved in SplineTextGeneration132.txt

n = length(RealminimalModel.TIME);
file = fopen('SplineTextGeneration132.txt','w');
fprintf(file,'number of knots: ' + string(n+2) + '\n\n');

fprintf(file,'TIME:\n[' + string(RealminimalModel.TIME(1)) + ', ');
for i = 1:n
    fprintf(file,string(RealminimalModel.TIME(i)) + ', ');
end
fprintf(file,string(RealminimalModel.TIME(n)) + ']\n\n');

fprintf(file,'INS:\n[' + string(RealminimalModel.INS(1)) + ', ');
for i = 1:n
    fprintf(file,string(RealminimalModel.INS(i)) + ', ');
end
fprintf(file,string(RealminimalModel.INS(n)) + ']\n\n');

