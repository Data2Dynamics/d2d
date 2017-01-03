% Progress bar for cluster using parfor_progress function
%
%   arShowProgressParFor(N) 
%       initializes the progress monitor for a set of N upcoming calculations.
%
%   arShowProgressParFor(j, n, startTime, custom_text) 
%       updates the progress inside your parfor loop and displays an updated.
%
%   arShowProgressParFor(0) deletes parfor_progress.txt and finalizes progress bar.

function arShowProgressParFor(j, n, startTime, custom_text)

thisworker = getCurrentWorker; % Worker object
if(~isempty(thisworker) && isfield(thisworker, 'Name'))
    worker_name = thisworker.Name;
else
    worker_name = 'local';
end

if(~exist('custom_text', 'var') || isempty(custom_text))
    custom_text = '';
end

if(nargin==1)
    pct = parfor_progress(j); %#ok<NASGU>
else 
    pct = parfor_progress/100;
    jpct = pct*n;
    timeelapsed = etime(clock, startTime);
    timeleft = (n-jpct) * timeelapsed/(jpct-1); % mean over all points
    
    fprintf('%i/%i [%s -> %s] fit #%i (%s, %s): %s\n', round(n*pct), n, ...
        secToHMS(timeleft), secToHMS(timeelapsed), j, ...
        worker_name, datestr(now, 0), custom_text);
end