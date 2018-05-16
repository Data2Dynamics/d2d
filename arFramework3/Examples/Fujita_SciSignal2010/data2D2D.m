clear all;
close all;
clc;

%% Data
% Experimental data is read out and written to an AMICI-data object which 
% is used for the ODE integration
D = getData_Akt_pathway();
obs_names = {'pEGFR_tot','pAkt_tot','pS6_tot'}; % define observable names
header = {'Time','EGF_conc_step','EGF_conc_impulse','EGF_conc_ramp','pulse_time','ramp_time','pEGFR_tot','pEGFR_tot_std','pAkt_tot','pAkt_tot_std','pS6_tot','pS6_tot_std'}; % define header for data files

%% Create D2D data format
for n = 1:length(D)
    tmp = [D(n).Y D(n).Sigma_Y];
    cell2csv(['./Data/experimentaldata' num2str(n) '.csv'],[header ; num2cell([D(n).t transpose(D(n).condition).*ones(length(D(n).t),length(D(n).condition)) tmp(:,[1 4 2 5 3 6])])])
    
    fileID = fopen(['./Data/experimentaldata' num2str(n) '.def'],'w');
    fprintf(fileID,['DESCRIPTION\n"Experiment ' num2str(n) '"\n\n']);
    fprintf(fileID,['PREDICTOR\nt\tT\t"min"\t"time"\t0\t' num2str(D(n).t(end)) '\n\n']);
    fprintf(fileID,'INPUTS\n\n');
    fprintf(fileID,'OBSERVABLES\n');
    for obs = 1:length(obs_names)
        fprintf(fileID,[ obs_names{obs} '\tC\tau\tconc.\t0\t0\t"' obs_names{obs} '_abs"\n' ]);
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'ERRORS\n');
    for obs = 1:length(obs_names)
        fprintf(fileID,[ obs_names{obs} '\t"sd_' obs_names{obs} '"\n' ]);
    end
    fprintf(fileID,'\n');
    fprintf(fileID,'CONDITIONS\n\n');
    fclose(fileID);
end