function int=setSeedPopFit(int, seedData, seedDataFit)
%   int=setSeedPopulation(int, seedData)
% Set the seed data for optimization. A 2-D array of doubles is expected
% which will be converted to the appropriate data type if necessary.
% The seed data must fit the problem dimension, while the population size
% will be adapted if necessary. If the seed data is an empty array, the
% seed population will be initialized randomly as by default.
%    int: the interface instance
%    seedData: 2-D array of dimension popSize x problemDim or [] to reset
%       and use random initial population

if length(seedData)>1
    %size(int.range,1)
    %int.dataType
    %int.dim
    if size(seedData,2)~=int.dim
        error(['Mismatching dimension: seed data should be of size popSize x problemDim. Current problemDim is ' num2str(int.dim)]);
    end
    if size(seedData,1)~=size(seedDataFit,1)
        error('Mismatching dimension of fitness array: expecting equal no. lines in seedData and seedDataFit!');
    end
    int.seedPop=seedData;
    int.seedPopFit=seedDataFit;
else
    int.seedPop=[]; % reset it
    int.seedPopFit=[];
end
