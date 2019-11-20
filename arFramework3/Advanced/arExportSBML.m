% arExportSBML(FileOption)
% 
% Exports current model to SBML
% Either as single SBML file (default) or with
% one SBML file for each condition
%
% FileOption:       'single' (default)
%                   'multi' (deprecated!)

function arExportSBML(varargin)
    global ar
    
    if (isempty(varargin))
        FileOption = 1;
    else
        if(varargin{1}=='single')
            FileOption=1;
        elseif(varargin{0}=='multi')
            FileOption=2;
        else
           error('Unclear function input. Please input single or multi!') 
        end
    end

    %Calculate chi2 value without Bessel correction
    ar.config.useFitErrorCorrection = false;
    if(contains(pwd,'Chen'))
        arCalcMerit
    else
        arFit
        arCalcMerit
    end

    if FileOption == 1
        for i = 1:length(ar.model)   
           arExportSBML_FullModel(i);       
        end
    elseif FileOption == 2
        for i = 1:length(ar.model)   
            for j = 1:length(ar.model(i).data)
        %           Export SBML files
                arExportSBML_datamodel(i,j,1);       
            end
        end
    end