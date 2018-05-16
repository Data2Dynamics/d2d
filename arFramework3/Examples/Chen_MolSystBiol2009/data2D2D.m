clear all;
close all;
clc;

%% Data
% Experimental data is read out and written to an AMICI-data object which 
% is used for the ODE integration

load('erbb_signaling_pnom.mat'); clear pnom;
D = getData_ErbB_signaling();
obs_names = {'ERB_B1_P_tot','ERK_PP','AKT_PP'}; % define observable names
header = {'Time','c1','c515','c285','c514','ERB_B1_P_tot','ERB_B1_P_tot_std','ERK_PP','ERK_PP_std','AKT_PP','AKT_PP_std'}; % define header for data files

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