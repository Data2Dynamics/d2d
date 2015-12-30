% Run a test for arShowProgressParFor
% either run script directly or send to cluster by:
%
% create cluster object, e.g.: 
%   c = parcluster('local');
%
% send batch job:
%   job = batch(c, 'arShowProgressParForTest', 'Pool', 3)

n = 100;

tic;

startTime = clock;
arShowProgressParFor(n);

parfor j=1:n
    pause(rand); % Replace with real code
    
    arShowProgressParFor(j, n, startTime, 'testing...');
end

arShowProgressParFor(0);

toc;