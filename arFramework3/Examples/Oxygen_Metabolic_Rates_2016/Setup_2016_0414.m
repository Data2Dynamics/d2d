%==========================================================================
% Setup file for PL analysis of 17O MRI data
%
%   updated 15.04.2016
%   
%==========================================================================
close all
clear

%==========================================================================
% chose experiment data 
%==========================================================================
N_exper = 1;    % 1 - res = 10 mm experiment, 2 - res = 8 mm experiment
BrainTissue = 1;    % 1 - WM, 2 - GM, 3 - CSF

%==========================================================================
% chose the model
% 1 - old model constant alpha and beta
% 2 - new model w prior area info
% 3 - new model w/o prior area info
% Note: alpha_dods can be used in N_model = 3 as prior. If no prior is
% assumed it should be commented out (line ~ 67 and 120)
%==========================================================================
N_model = 3;    % 

%==========================================================================
% load model and data and run fit and use known information
% first experiment (res = 10mm)
%==========================================================================
arInit
if N_exper == 1
    if N_model == 1
        arLoadModel('CMRO_Hoff_Exp1');
    elseif N_model == 2
        arLoadModel('CMRO_Exp1_V03');
        arLoadData('Input_Area_Exp2_1304'); % 1.33
    elseif N_model == 3
        arLoadModel('CMRO_Exp1_V03_wo_Prior');
    else
        error ('no proper model is chosen')
    end
    if BrainTissue == 1
        arLoadData('CMRO_data_WM_10mm_woH');
    elseif BrainTissue == 2
        arLoadData('CMRO_data_GM_10mm_woH');
    elseif BrainTissue == 3
        arLoadData('CMRO_data_CSF_10mm_woH');
    else
        error ('no 17O MRI data is chosen')
    end
    
    arCompileAll;
    arFindInputs;
    arSetPars('rho',0.75,2,0,-5,3);
    arSetPars('time_delay',10.5,2,0,9,12);  % correct for 1st measurement
    
    if N_model == 1
    	alpha_error = 10; % uncertanty in %
        arSetPars('alpha_input',0.27,1,0,0,1,1,0.27,0.27*alpha_error*0.01);
    end
    if N_model == 2
    	alpha_error = 10; % uncertanty in %
    	arSetPars('sd_Input_Flaeche',log10(1.33*alpha_error*0.01),2,1,-5,3);
    end
    if N_model == 3
    	alpha_error = 10; % uncertanty in %
    	arSetPars('alpha_input',0.27,1,0,0,1,1,0.27,0.27*alpha_error*0.01);
    end
    arFit
    arPlot
    arPrint
    
%     arPLEInit
%     ple(1:9)

end


%==========================================================================
% load model and data and run fit and use known information
% second experiment (res = 8mm)
%==========================================================================
if N_exper == 2
    arInit
    if N_model == 1
        arLoadModel('CMRO_Hoff_Exp2');
    elseif N_model == 2
        arLoadModel('CMRO_Exp2_V03');
        arLoadData('Input_Area_Exp2_1304'); % 1.29
    elseif N_model == 3
        arLoadModel('CMRO_Exp2_V03_wo_Prior');
    else
        error ('no proper model is chosen')
    end
    if BrainTissue == 1
        arLoadData('CMRO_data_WM_8mm_woH');
    elseif BrainTissue == 2
        arLoadData('CMRO_data_GM_8mm_woH');
    elseif BrainTissue == 3
        arLoadData('CMRO_data_CSF_8mm_woH');
    else
        error ('no 17O MRI data is chosen')
    end
    
    arCompileAll;
    arFindInputs;
    arSetPars('rho',0.75,2,0,-5,3);
    arSetPars('time_delay',9.17,2,0,9,12);  % correct for 2nd measurement
    if N_model == 1
    	alpha_error = 10; % uncertanty in %
    	arSetPars('alpha_input',0.31,1,0,0,1,1,0.31,0.31*alpha_error*0.01);
        % arSetPars('beta_input',log10(0.4765*0.5),2,1,-3,5);
    end
	if N_model == 2
    	alpha_error = 10; % uncertanty in %
    	arSetPars('sd_Input_Flaeche',log10(1.29*alpha_error*0.01),2,1,-5,3);
    end
    if N_model == 3
    	alpha_error = 10; % uncertanty in %
        arSetPars('alpha_input',0.31,1,0,0,1,1,0.31,0.31*alpha_error*0.01);
    end

    arFit
    arPlot
    arPrint
    
%     arPLEInit
%     ple(1:9)

end



