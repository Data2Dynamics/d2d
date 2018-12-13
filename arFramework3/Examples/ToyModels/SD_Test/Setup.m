arInit

arLoadModel('SD_model');
arLoadData('ABC_data_BCobs'); %Data with equidistant observation
%of state B and C for t=0,10,..100
%arLoadData('ABC_data_B_sparseObs'); %Data with sparse observation of
                                 %state B
arCompileAll();
