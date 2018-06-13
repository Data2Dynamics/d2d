%% loading model and a wild-type data-set where everything is observed
Setup

%% Loading of the original parameters and simulation of data:
arSimuData
arPlot

%% Simulating all single pertubations:
indper = [find(~cellfun(@isempty,regexp(ar.pLabel,'_ko$'))),find(~cellfun(@isempty,regexp(ar.pLabel,'_kd$'))),find(~cellfun(@isempty,regexp(ar.pLabel,'_ic$')))];
ps = ones(length(indper),1)*ar.p;
for i=1:length(indper)
    ps(i,indper(i)) = 1;
end
arPlotMulti(ps)



