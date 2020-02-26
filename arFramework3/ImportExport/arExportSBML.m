% arExportSBML([name], [ModelForEachCondition])
% 
% Exports current model to SBML
% Either as single SBML file (default) or with
% one SBML file for each condition
%
% name              string of name for model export
%                   default: name of current folder
%                   
% ModelForEachCondition:   
%                   '0'     Export model in single xml file('default')
%                   '1'     Export model for each condition 
%                   

function arExportSBML(name,ModelForEachCondition)
    global ar
      
    if ~exist('name','var') || isempty(name)
        directory = split(pwd,filesep);
        name = directory{end};
    end 
    if ~exist('ModelForEachCondition','var') || ...
            (ModelForEachCondition~=0 && ModelForEachCondition~=1)
        ModelForEachCondition = 0;
    end

    %Calculate chi2 value without Bessel correction
    ar.config.useFitErrorCorrection = false;
    disp('Deactivating Bessel correction...')
    arCalcMerit;

    if ModelForEachCondition == 0 % single file output
        for i = 1:length(ar.model)   
           arExportSBML_FullModel(i,name);       
        end
    elseif ModelForEachCondition == 1 % multi file output
        for i = 1:length(ar.model)   
            for j = 1:length(ar.model(i).data)
        %           Export SBML files
                arExportSBML_singlecondition(i,j,name);       
            end
        end
    end