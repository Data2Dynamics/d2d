% arExportSBML(FileOptionString, name)
% 
% Exports current model to SBML
% Either as single SBML file (default) or with
% one SBML file for each condition
%
% FileOptionString:       
%                   'multi' (default!)
%                   'single' ('not yet finished due to structural problems in d2d')

function arExportSBML(FileOptionString, name)
    global ar
    
    if ~exist('FileOption') || isempty(FileOption)
        FileOptionString = 'multi';
    end
    if ~exist('name') || isempty(name)
        name = 'name';
    end
    
    if(strcmp(FileOptionString,'single'))
        FileOption=1;
    elseif(strcmp(FileOptionString,'multi'))
        FileOption=2;
    else
        error('Unclear function input. Please input single or multi!')
    end

    %Calculate chi2 value without Bessel correction
    ar.config.useFitErrorCorrection = false;
    disp('Deactivating Bessel correction...')
    arCalcMerit;

    if FileOption == 1 % single file output
        for i = 1:length(ar.model)   
           arExportSBML_FullModel(i,[],name);       
        end
    elseif FileOption == 2 % multi file output
        for i = 1:length(ar.model)   
            for j = 1:length(ar.model(i).data)
        %           Export SBML files
                arExportSBML_singlecondition(i,j,1,name);       
            end
        end
    end