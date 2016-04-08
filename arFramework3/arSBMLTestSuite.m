
% Import and simulate all SBML Test Suite models (http://sbml.org/Facilities/Database) and
% compare results with reference simulation.

d2d_dir = fileparts(which('arInit.m'));

if(~exist([d2d_dir '/sbml-test-cases-2014-10-22'],'dir'))
    fprintf(1, '\nDownloading SBML test cases...\n');
    websave([d2d_dir '/sbml-test-cases-2014-10-22.zip'],'http://jaist.dl.sourceforge.net/project/sbml/test-suite/3.1.1/cases-archive/sbml-test-cases-2014-10-22.zip');
    fprintf(1, '\nUnizpping SBML test cases...\n');
    unzip([d2d_dir '/sbml-test-cases-2014-10-22.zip'],[d2d_dir '/sbml-test-cases-2014-10-22']);
end

if ~exist('TestSuite','dir')
    mkdir('TestSuite')
end

cd('TestSuite')

modeldir = [d2d_dir '/sbml-test-cases-2014-10-22/cases/semantic/'];
dirname = dir(modeldir);
dirname = dirname([dirname.isdir]);
dirname = dirname(3:end);

results_table = array2table({'max_deviation','error_msg'});
results_table.Properties.VariableNames = results_table{1,:};
results_table(1,:) = [];

for jm = 1:length(dirname)

    modelnum = dirname(jm).name;
    settings = readtable([modeldir modelnum '/' modelnum '-settings.txt'],'Delimiter',':',...
        'ReadVariableNames',false,'ReadRowNames',true);

    xml_name = [modeldir modelnum '/' modelnum '-sbml-l3v1.xml'];
    copyfile(xml_name, [modelnum '-sbml-l3v1.xml'])


    modelname = [modelnum '-sbml-l3v1'];
    tmax = str2num(char(settings{'duration',:}));
    nsteps = str2num(char(settings{'steps',:}));
    
    atol = str2num(char(settings{'absolute',:}));
    rtol = str2num(char(settings{'relative',:}));
    
    arInit
    % ar.config.nFinePoints = nsteps;
    ar.config.nFinePoints = 500;
    % ar.config.atol = atol;
    % ar.config.rtol = rtol;
    ar.config.checkForNegFluxes = false;
    try
        [~, modelname] = arImportSBML(modelname,tmax);
        arLoadModel(modelname)
        arCompileAll

        arSimu(0,1,1);

        results = readtable([modeldir dirname(jm).name '/' modelnum '-results.csv'],'ReadVariableNames',true,'ReadRowNames',false);
        varnames = results.Properties.VariableNames(2:end);
        
        % find variables and inputs (boundarySpecies)
        xNames = intersect(varnames,ar.model.xNames);
        uNames = intersect(varnames,ar.model.u);
        
        % interpolate time
        d2d_res = interp1(ar.model.condition.tFine,[ar.model.condition.xFineSimu ar.model.condition.uFineSimu],results.time);
        
        deviation = abs(results{:,[xNames uNames]}-d2d_res);
        results_table{modelnum,:} = {max(deviation(:)), []};
        
        if all(deviation(:)<1e-5))
            fprintf(1, '\n\nModel %s sucessfully imported! Maximal deviation: %2.2g\n\n', modelnum, max(deviation(:)));
        else
            warning('Maximal deviation >1e-5');
            clf
            plot(results.time,results{:,2:end},'-')
            hold on
            plot(results.time,d2d_res,'o--')
            hold off
            title(modelnum)
            legend([varnames varnames])
            pause    
        end
    catch err
        results_table{modelnum,:} = {[], err.message};
    end
end

% print result_table
results_table %#ok<NOPTS>
