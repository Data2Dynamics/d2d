clear all;
close all;
clc;

fprintf('Loading AMICI model ... ')
model = erbb_signaling_pesto_syms(); % execute AMICI model function
obs_names = {'ERB_B1_P_tot','ERK_PP','AKT_PP'}; % define observable names
fprintf('done!\n')

%% Create model file for D2D from AMICI
fprintf('Writing D2D model file ... ')
fileID = fopen('./Models/erbb_signaling.def','w');

fprintf(fileID,['DESCRIPTION\n"Model of the ErbB signaling pathways"\n']);
fprintf(fileID,['"This model is taken from the paper: Input-output behavior of ErbB signaling pathways as revealed by a mass action model trained against dynamic data, by Chen et al., in 2009, PubMed, vol.5, no.239"\n']);
fprintf(fileID,['"Available at https://www.ncbi.nlm.nih.gov/pubmed/19156131"\n\n']);

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

fprintf(fileID,['OBSERVABLES\ndummi\tC\tau\tconc.\t0\t0\t"1"\n\n']);

fprintf(fileID,['ERRORS\ndummi\t"1"\n\n']);

fprintf(fileID,'CONDITIONS\n');
for n = 1:length(model.sym.x)
    fprintf(fileID,['init_' char(model.sym.x(n)) '\t"' num2str(model.sym.x0(n)) '"\n' ]);
end
fclose(fileID);

fprintf('done!\n')