% waitbar with time estimation

function arWaitbar(j, n, text)

global arWaitbarGlobal;

% disable waitbar window and use command line output instead
arWaitbar.showWindow = 1;

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
    if(isfield(arWaitbarGlobal,'h'))
        try %#ok<TRYNC>
            close(arWaitbarGlobal.h);
        end
    end
    
    arWaitbarGlobal.tic = clock;
    arWaitbarGlobal.ticlast = clock;
    arWaitbarGlobal.ticmove = clock;
    monitorPositions = get(0,'MonitorPositions');
    arWaitbarGlobal.hasMonitor = sum(sum(repmat([0 0 1 1],size(monitorPositions,1),1) == monitorPositions,2) ~= 4) > 0;
    
    arWaitbarGlobal.timeperstep = [];
    
    arWaitbarGlobal.tcount = 1;
	
    if(arWaitbarGlobal.hasMonitor && arWaitbar.showWindow)
        if(exist('text','var') && ~isempty(text))
            arWaitbarGlobal.h = waitbar(0, text);
        else
            arWaitbarGlobal.h = waitbar(0, 'Please wait...');
        end
    end
elseif(j<0)
    if(arWaitbarGlobal.hasMonitor && arWaitbar.showWindow)
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

            timeleft = (n-j) * timeelapsed/(j-1); % mean over all points
            
            funtext = '';
            if(timeleft > 60*5)
                funtext = '   (get a coffee !)';
            end
            if(timeleft > 60*30)
                funtext = '   (go for lunch !)';
            end
            if(timeleft > 60*60*2)
                funtext = '   (go home !)';
            end
            if(timeleft > 60*60*24*7)
                funtext = '   (come back next week !)';
            end
            if(timeleft > 60*60*24*7*30)
                funtext = '   (this must be a joke !)';
            end

            if(exist('text','var') && ~isempty(text))
                if(arWaitbarGlobal.hasMonitor && arWaitbar.showWindow)
                    arWaitbarGlobal.h = waitbar(j/n, arWaitbarGlobal.h, ...
                        sprintf('%s\n%i/%i   %2i%%   %s -> %s%s', ...
                        text, j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed), funtext));
                else
                    fprintf('%s %i/%i   %2i%%   %s -> %s%s', ...
                        text, j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed), funtext); 
                    drawnow;
                end
            else
                if(arWaitbarGlobal.hasMonitor && arWaitbar.showWindow)
                    arWaitbarGlobal.h = waitbar(j/n, arWaitbarGlobal.h, ...
                        sprintf('%i/%i   %2i%%   %s -> %s%s', ...
                        j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed), funtext));
                else
                    fprintf('%i/%i   %2i%%   %s -> %s%s\n', ...
                        j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed), funtext); 
                    drawnow;
                end
            end
            
            arWaitbarGlobal.ticlast = clock;
        end
    end
end
