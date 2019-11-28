% arExportSBML(FileOption)
% 
% Exports current model to SBML
% Either as single SBML file (default) or with
% one SBML file for each condition
%
% FileOption:       
%                   'multi' (default!)
%                   'single' ('not yet finished due to structural problems in d2d')

function arExportSBML(varargin)
    global ar
    
    if (isempty(varargin))
        FileOption = 2;
    else
        if(strcmp(varargin{1},'single'))
            FileOption=1;
        elseif(strcmp(varargin{1},'multi'))
            FileOption=2;
        else
           error('Unclear function input. Please input single or multi!') 
        end
    end

    %Calculate chi2 value without Bessel correction
    ar.config.useFitErrorCorrection = false;
    arCalcMerit;

    if FileOption == 1 % single file output
        for i = 1:length(ar.model)   
           arExportSBML_FullModel(i);       
        end
    elseif FileOption == 2 % multi file output
        for i = 1:length(ar.model)   
            for j = 1:length(ar.model(i).data)
        %           Export SBML files
                arExportSBML_singlecondition(i,j,1);       
            end
        end
    end