
function results_table = arSBMLTestSuite(model_id, sbml_ver, save_output, stoponwarning)
% Import and simulate all SBML Test Suite models (http://sbml.org/Facilities/Database) and
% compare results with reference simulation.
% 
% results_table = arSBMLTestSuite(model_id, sbml_ver)
% model_id: 1:1196
% sbml_ver: 'l3v1'

if(~exist('model_id','var') || isempty(model_id))
    model_id = [];
end
if(~exist('sbml_ver','var') || isempty(sbml_ver))
    sbml_ver = 'l3v1';
end
if(~exist('save_output','var') || isempty(save_output))
    save_output = false;
end
if(~exist('stoponwarning','var') || isempty(stoponwarning))
    stoponwarning = true;
end

sbml_ver = ['-sbml-' sbml_ver];

d2d_dir = fileparts(which('arInit.m'));

if(~exist([d2d_dir '/sbml-test-cases-2014-10-22.zip'],'file'))
    fprintf(1, '\nDownloading SBML test cases...\n');
    websave([d2d_dir '/sbml-test-cases-2014-10-22.zip'],'http://jaist.dl.sourceforge.net/project/sbml/test-suite/3.1.1/cases-archive/sbml-test-cases-2014-10-22.zip');
    if(~exist([d2d_dir '/sbml-test-cases-2014-10-22'],'dir'))
        fprintf(1, '\nUnizpping SBML test cases...\n');
        unzip([d2d_dir '/sbml-test-cases-2014-10-22.zip'],[d2d_dir '/sbml-test-cases-2014-10-22']);
    end
end

if ~exist('TestSuite','dir')
    mkdir('TestSuite')
    cd('TestSuite')
end

modeldir = [d2d_dir '/sbml-test-cases-2014-10-22/cases/semantic/'];
dirname = dir(modeldir);
dirname = dirname([dirname.isdir]);
dirname = dirname(3:end);

results_table = array2table({'max_deviation','error_msg'});
results_table.Properties.VariableNames = results_table{1,:};
results_table(1,:) = [];

if(isempty(model_id))
    model_id = 1:length(dirname);
end

for jm = model_id
    modelnum = dirname(jm).name;
    settings = readtable([modeldir modelnum '/' modelnum '-settings.txt'],'Delimiter',':',...
        'ReadVariableNames',false,'ReadRowNames',true);

    xml_name = [modeldir modelnum '/' modelnum sbml_ver '.xml'];
    copyfile(xml_name, [modelnum sbml_ver '.xml'])


    modelname = [modelnum sbml_ver];
    tmax = str2num(char(settings{'duration',:}));
    nsteps = str2num(char(settings{'steps',:}));
    
    % atol = str2num(char(settings{'absolute',:}));
    % rtol = str2num(char(settings{'relative',:}));
    
    arInit
    ar.config.nFinePoints = nsteps+1;
    ar.config.atol = 1e-9;%atol;
    ar.config.rtol = 1e-9;%rtol;
    ar.config.checkForNegFluxes = false;
    
    results = readtable([modeldir dirname(jm).name '/' modelnum '-results.csv'],'ReadVariableNames',true,'ReadRowNames',false);
    varnames = results.Properties.VariableNames(2:end);

    try
        [m, modelname] = arImportSBML(modelname,tmax);
        arLoadModel(modelname)
        arCompileAll
        arSimu(0,1,1);
        
        % find variables and inputs (boundarySpecies)
        xNames = intersect(varnames,ar.model.xNames);
        uNames = intersect(varnames,ar.model.u);
        
        % collect d2d results
%         if all([m.d2d.species(:).isSetInitialAmount]==1)
%             d2d_tmp = [ar.model.condition.xFineSimu ar.model.condition.uFineSimu]*0.3;
%         else
            d2d_tmp = [ar.model.condition.xFineSimu ar.model.condition.uFineSimu];
%         end
        
        d2d_results = results;
        d2d_results{:,[xNames uNames]} = d2d_tmp;
        if save_output
            writetable(d2d_results,['d2d_' modelnum '.csv'],'Delimiter',',');
        end
        
        deviation = abs(results{:,:}-d2d_results{:,:});
        results_table{modelnum,:} = {max(deviation(:)), []};
        
        if all(deviation(:)<1e-5)
            fprintf(1, '\n\nModel %s sucessfully imported! Maximal deviation: %2.2g\n\n', modelnum, max(deviation(:)));
        elseif stoponwarning
            warning('Maximal deviation >1e-5');
            
            clf;
            plot(results.time,results{:,2:end},'-')
            hold on
            plot(results.time,d2d_results{:,2:end},'o--')
            hold off
            title(modelnum)
            legend([varnames varnames])
            
            fprintf(1,'\nPress any key to continue\n\n');
            pause;
            
            close gcf;
        end
    catch err
        results_table{modelnum,:} = {[], err.message};
        delete('*.def')
    end
    
end

