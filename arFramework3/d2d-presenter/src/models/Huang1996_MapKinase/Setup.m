%% Create *.def
arImportSBML('BIOMD0000000009_d2d',100)

%% Load models & data
arInit
ar.config.checkForNegFluxes = false
arLoadModel('BIOMD0000000009_d2d');
arLoadData('BIOMD0000000009_d2d_data');
arCompileAll;

%% Check
ar.config.atol=1e-10;
ar.config.rtol=1e-10;
%arQplot('x')
%arPlot

%arCompareWithBiobaseSimulation('SIMU1447942223190.dat');

%%
e1vals = logspace(-6,-1,101);
out = NaN(length(e1vals),3);
ix = [4,7,10];
for i=1:length(e1vals)
    ar.p(1) = e1vals(i);
%     arSimu(false,true,true);
    arSimu(false, true);
    arSimu(false, false);
    out(i,:) = ar.model.condition.xFineSimu(end,ix);
end

%out = out./(ones(size(out,1),1)*max(out,[],1));
%semilogx(e1vals,out)
%legend(strrep(ar.model.x(ix),'_',' '),'Location','SouthEast');