function int=setSeedPopulation(int, seedData)
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
    %size(seedData)
    %size(int.range,1)
    %int.dataType
    %int.dim
    if size(seedData,2)~=int.dim
        error(['Mismatching dimension: seed data should be of size popSize x problemDim. Current problemDim is ' num2str(int.dim)]);
    end
    
    for i=1:size(seedData,1)
        if (isempty(int.args))
            fit(i,:) = feval(int.f, seedData(i,:));
        else
            fit(i,:) = feval(int.f, seedData(i,:), int.args);
        end
    end
    int=setSeedPopFit(int,seedData,fit);
else
    int.seedPop=[]; % reset it
    int.seedPopFit=[];
end
