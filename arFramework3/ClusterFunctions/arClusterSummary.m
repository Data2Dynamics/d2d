function arClusterSummary(c)

fprintf('%s: %i workers, %i busy, %i idle\n', ...
    c.Profile, ...
    c.NumWorkers, c.NumBusyWorkers, ...
    c.NumIdleWorkers);

[pending, queued, running, finished] = findJob(c);
subprint(pending, 'pending  ', 'SubmitTime', 1);
subprint(queued, 'queued   ', 'SubmitTime', 2);
subprint(running, 'running  ', 'StartTime', 3);
subprint(finished, 'finished ', 'FinishTime', 4);



function subprint(jobs, str, timefield, tid)
fprintf('%s: %2i jobs\n', str, length(jobs));
for j=1:length(jobs)
    timestr = jobs(j).(timefield);
    if(tid==1 || tid==2 || tid==3)
        fid = 1;
    elseif(tid==4)
        if(isempty(jobs(j).Tasks(1).ErrorIdentifier))
            fid = 1;
        else
            fid = 2;
        end
    end
    fprintf(fid, '#%i ID %i: %2i tasks (%2i-%2i), %s %s, %s\n', j, jobs(j).ID, ...
        length(jobs(j).Tasks), jobs(j).NumWorkersRange, timefield, timestr, jobs(j).Username);
end
