%Testing x,y,z Fine Simus and sx, sy, sz Exp Simus for all Benchmark models
%published at XXX
arInit;

list_examples = {'Becker_Science2010';'Bachmann_MSB2011';'Beer_MolBiosyst2014';'Boehm_JProteomeRes2014';'Bruno_Carotines_JExpBio2016';'TGFb_ComplexModel_WithGenes_Reduced'; ...
    'Merkle_JAK2STAT5_PCB2016';'Raia_CancerResearch2011';'Reelin_PONE2017';'Schwen_InsulinMouseHepatocytes_PlosOne2014';'Swameye_PNAS2003'; ...
'Brannmark_JBC2010';'Crauste_ImmuneCells_CellSystems2017';'Weber_BMC2015';'Isensee_JCB2018';'Zheng_PNAS2012';'Fiedler_BMC2016';'Sobotta_Frontiers2017';'Fujita_SciSignal2010'};%;'Chen_MSB2009'

setups_examples = {'Setup','Setup','Setup_IndiBac','Setup_FullModel_Boehm2014','Setup','Setup','SetupFinal','Setup','Setup_final','Setup','Setup','Setup','Setup','Setup','Setup','Setup_Zheng','Setup','Setup_Core','Setup'};%,'Setup'

load('Test_Benchmark_Simulations.mat');
list_fields = {'yFineSimu','syExpSimu','xFineSimu','sxExpSimu','zFineSimu','szExpSimu'};

diff_error = false;
for iEx = 1:length(list_examples)
    if isunix
        command = sprintf('find %s -type d -name "%s"',example_folder,list_examples{iEx});
    else
        command = sprintf('dir "%s" /AD /b /s | findstr /r "\\\\%s$"',example_folder,list_examples{iEx});
    end
    [~,dir] = system(command);
    dir = strread(dir, '%s', 'delimiter', sprintf('\n'));
    cd(dir{1});
   
   did_load = arLoadLatest('compiled');
    
   if(~did_load)
       sprintf('Recompiling %s \n',list_examples{iEx})
       eval(setups_examples{iEx})
       arSave('compiled')
   end
   fprintf('Loaded model %s \n', list_examples{iEx})
   close all
   %Initialize same nFinePoints
   ar.config.nFinePoints = 100;
   arLink
   
   arSimu(1,0,1)
   arSimu(0,1,1)
   
   
   for iFields = 1:length(list_fields)
      if(list_fields{iFields}(1) == 'y' || list_fields{iFields}(2) == 'y')
          diff_tmp = Simulation_struct.(list_examples{iEx}).(list_fields{iFields}) - ar.model(1).data(1).(list_fields{iFields});
      else
          diff_tmp = Simulation_struct.(list_examples{iEx}).(list_fields{iFields}) - ar.model(1).condition(1).(list_fields{iFields});
      end
      if(diff_tmp > 1.e-4)
          diff_error = true;
      end
   end
   if(~diff_error)
       fprintf( 'PASSED %s \n', list_examples{iEx} );
   else
       error( 'ERROR BETWEEN SIM AND STORED VALUE TOO LARGE' );
    end
end 


