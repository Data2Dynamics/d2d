% arSimuMeasurements
% 
% arSimuMeasurements(preserveZeros,varargin)
% 
% Simulates data like arSimuData, but preserves availablity of data as in the
% original data files. In other words, if NaN (or 0) is in
% ar.model.data.yExpRaw, then this is NOT replaced by a random number.
% 
% Like for arSimu, the simulated data has the "right" std.dev, i.e. the
% noise level is given by the prespecified exp. error or by the error
% model.
% 
%   preserveZeros   default: 1 (zeros are preserved)
%                   alternative: 0
% 
%   varargin        arguments passed to arSimu

function arSimuMeasurements(preserveZeros,varargin)
if ~exist('preserveZeros','var') || isempty(preserveZeros)
    preserveZeros = 1;
end
    
global ar

arSimuData(varargin{:}) % simuliert auch für yExp=0 und yExp=NaN

for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        if preserveZeros
            ar.model(m).data(d).yExp(ar.model(m).data(d).yExpRaw==0) = 0;
        end
        ar.model(m).data(d).yExp(isnan(ar.model(m).data(d).yExpRaw)) = NaN;
    end
end

