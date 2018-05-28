arInit;

arLoadModel('model_template');

arLoadData('data_template1');
arLoadData('data_template2');

arCompileAll;

ar.config.useNewPlots=1;
ar.config.fiterrors=1;

ar.p(1) = 0;
ar.p(2) = log10(.2);
ar.p(3) = log10(.5);
ar.p(4) = log10(.5);
ar.p(5) = log10(2);
ar.p(6:7) = 0;
ar.p(8) = log10(.2);

