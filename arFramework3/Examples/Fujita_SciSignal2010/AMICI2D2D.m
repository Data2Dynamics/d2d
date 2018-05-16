clear all;
close all;
clc;

fprintf('Loading AMICI model ... ')
model = AktPathway_syms(); % execute AMICI model function
obs_names = {'pEGFR_tot','pAkt_tot','pS6_tot'}; % define observable names
fprintf('done!\n')

%% Create model file for D2D from AMICI
fprintf('Writing D2D model file ... ')
fileID = fopen('./Models/AktPathwayFujita.def','w');

fprintf(fileID,['DESCRIPTION\n"Model of the PI3K-EGFR signaling pathway"\n']);
fprintf(fileID,['"This model is taken from the paper "Decoupling of Receptor and Downstream Signals in the Akt Pathway by Its Low-Pass Filter Characteristics" , by Fujita et. al., in 2010, ScienceSignaling, Vol.3 Issue 132 ra56"\n']);
fprintf(fileID,['"https://www.ncbi.nlm.nih.gov/pubmed/20664065"\n\n']);

fprintf(fileID,['PREDICTOR\nt\tT\t"min"\t"time"\t0\t240\n\n']);

fprintf(fileID,['COMPARTMENTS\ncyt\tV\t"pl"\t"vol."\t1\n\n']);

fprintf(fileID,'STATES\n');
for n = 1:length(model.sym.x)
    fprintf(fileID,[char(model.sym.x(n)) '\tC\tmumol/l\tconc.\tcyt\t1\t"' char(model.sym.x(n)) '"\n' ]);
end
fprintf(fileID,'\n');

fprintf(fileID,'INPUTS\n');
for n = 1:length(model.sym.k)
    fprintf(fileID,[char(model.sym.k(n)) '\tC\tau\tconc.\t"' char(model.sym.k(n)) '"\n' ]);
end
fprintf(fileID,'\n');

fprintf(fileID,'ODES\n');
for n = 1:length(model.sym.xdot)
    fprintf(fileID,['"' char(model.sym.xdot(n)) '"\n' ]);
end
fprintf(fileID,'\n');

fprintf(fileID,'DERIVED\n');
for n = 1:length(model.sym.y)
    fprintf(fileID,[obs_names{n} '_abs\tC\tmumol/l\tconc.\t"' char(model.sym.y(n)) '"\n' ]);
end
fprintf(fileID,'\n');

% fprintf(fileID,['OBSERVABLES\ndummi\tC\tau\tconc.\t0\t0\t"1"\n\n']);
fprintf(fileID,['OBSERVABLES\n\n']);

% fprintf(fileID,['ERRORS\ndummi\t"1"\n\n']);
fprintf(fileID,['ERRORS\n\n']);

fprintf(fileID,'CONDITIONS\n');
for n = 1:length(model.sym.x)
    if isnumeric(model.sym.x0(n))
        fprintf(fileID,['init_' char(model.sym.x(n)) '\t"' num2str(model.sym.x0(n)) '"\n' ]);
    else
        fprintf(fileID,['init_' char(model.sym.x(n)) '\t"' char(model.sym.x0(n)) '"\n' ]);
    end
end
fclose(fileID);

fprintf('done!\n')