% waitbar with time estimation

function arWaitbar(j, n, text)

global arWaitbarGlobal;

% Ntimes = 20;

if(length(j) > 1)
    if(length(j) == length(n))
        jsum = j(end);
        for jj = (length(n)-1):-1:1
            jsum = jsum + (j(jj)-1)*prod(n((jj+1):end));
        end
        j = jsum;
        n = prod(n);
    else
        error('length(j) != length(n)');
    end
end

if(j==0 && nargin==1)
    arWaitbarGlobal.tic = clock;
    arWaitbarGlobal.ticlast = clock;
    arWaitbarGlobal.ticmove = clock;
    arWaitbarGlobal.hasMonitor = sum(get(0,'MonitorPositions')==[0 0 1 1]) ~= 4;
    
    arWaitbarGlobal.timeperstep = [];
    
    arWaitbarGlobal.tcount = 1;
	
    if(arWaitbarGlobal.hasMonitor)
        if(exist('text','var') && ~isempty(text))
            arWaitbarGlobal.h = waitbar(0, text);
        else
            arWaitbarGlobal.h = waitbar(0, 'Please wait...');
        end
    end
elseif(j<0)
    if(arWaitbarGlobal.hasMonitor)
        try %#ok<TRYNC>
            close(arWaitbarGlobal.h);
        end
    end
else
    try %#ok<TRYNC>
        
        if(j==1)
            arWaitbarGlobal.timeperstep = nan(1,n);
        else
            arWaitbarGlobal.timeperstep(j) = etime(clock, arWaitbarGlobal.ticmove);
        end
        arWaitbarGlobal.ticmove = clock;
        
        if(etime(clock, arWaitbarGlobal.ticlast)>arWaitbarGlobal.tcount)
            timeelapsed = etime(clock, arWaitbarGlobal.tic);

            timeleft = (n-j) * (timeelapsed/j); % mean over all points

%             js = find(~isnan(arWaitbarGlobal.timeperstep));
%             if(length(js)>Ntimes)
%                 js = js((end-Ntimes):end);
%             end
%             times = arWaitbarGlobal.timeperstep(js);
%             js = js(~isnan(times));
%             times = times(~isnan(times));
            
%             timeleft = (n-j) * median(times); % mean over last points
            
%             mjs = mean(js);
%             mtimes = mean(times);
%             js = [js js];
%             times = [times mtimes*ones(size(times))];
%             b = sum((times-mtimes).^2)/sum((js-mjs).*(times-mtimes));
%             a = mtimes - b*mjs;
%             predtime = a + b.*((j+1):n);
%             timeleft = sum(predtime(predtime>0));
%             if(timeleft<0)
%                 timeleft = 0;
%             end
            
            if(exist('text','var') && ~isempty(text))
                if(arWaitbarGlobal.hasMonitor)
                    arWaitbarGlobal.h = waitbar(j/n, arWaitbarGlobal.h, ...
                        sprintf('%s\n%i/%i | %2i%% | %s -> %s', ...
                        text, j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed)));
                else
                    fprintf('%s %i/%i | %2i%% | %s -> %s\n', ...
                        text, j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed)); 
                    drawnow;
                end
            else
                if(arWaitbarGlobal.hasMonitor)
                    arWaitbarGlobal.h = waitbar(j/n, arWaitbarGlobal.h, ...
                        sprintf('%i/%i | %2i%% | %s -> %s', ...
                        j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed)));
                else
                    fprintf('%i/%i | %2i%% | %s -> %s\n', ...
                        j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed)); 
                    drawnow;
                end
            end
            
%             arWaitbarGlobal.timeperstep(j) = nan;
            arWaitbarGlobal.ticlast = clock;
        end
    end
end