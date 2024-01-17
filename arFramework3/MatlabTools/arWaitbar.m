% arWaitbar(j, n, text)
% waitbar with time estimation

function arWaitbar(j, n, text)

persistent batchmode
if isempty(batchmode)
    ss = get(0, 'ScreenSize');
    if max(ss)<2
        batchmode = true;
    else
        batchmode = false;
    end
end

if batchmode 
    return % do not show Waitbar in batch mode
end

global arWaitbarGlobal;
global arOutputLevel;

global ar;		
% disable waitbar window and use command line output instead
if ( isfield( ar, 'config' ) && isfield( ar.config, 'noWaitBarWindow' ) && ar.config.noWaitBarWindow == 1 )
    arWaitbarGlobal.showWindow = 0;
else
    if ~isfield(arWaitbarGlobal,'showWindow')
        arWaitbarGlobal.showWindow = 1;
    end
end

if ( isfield( ar, 'config' ) && isfield( ar.config, 'noWaitBar' ) && ar.config.noWaitBar == 1 )		
    return;		
end

% suppress waitbar
if ( arOutputLevel < 2 )
    return;
end

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
	
    if(arWaitbarGlobal.hasMonitor && arWaitbarGlobal.showWindow)
        if(exist('text','var') && ~isempty(text))
            arWaitbarGlobal.h = waitbar(0, strrep(text, '_','\_'));
        else
            arWaitbarGlobal.h = waitbar(0, 'Please wait...');
        end
    end
elseif(j<0)
    if(arWaitbarGlobal.hasMonitor && arWaitbarGlobal.showWindow)
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
            if(timeleft > 60*10)
                funtext = '   (get a coffee !)';
            end
            if(timeleft > 60*60)
                funtext = '   (go for lunch !)';
            end
            if(timeleft > 60*60*3)
                funtext = '   (go home !)';
            end
            if(timeleft > 60*60*24*7)
                funtext = '   (come back next week !)';
            end
            if(timeleft > 60*60*24*7*30)
                funtext = '   (are you serious ?)';
            end

            if(exist('text','var') && ~isempty(text))
                if(arWaitbarGlobal.hasMonitor && arWaitbarGlobal.showWindow)
                    arWaitbarGlobal.h = waitbar(j/n, arWaitbarGlobal.h, ...
                        sprintf('%s\n%i/%i   %2i%%   %s -> %s%s', ...
                        strrep(text, '_','\_'), j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed), funtext));
                else
                    fprintf('%s %i/%i   %2i%%  Estimated: %s -> Elapsed: %s%s\n', ...
                        strrep(text, '_','\_'), j, n, round(j/n*100), secToHMS(timeleft), secToHMS(timeelapsed), funtext); 
                    drawnow;
                end
            else
                if(arWaitbarGlobal.hasMonitor && arWaitbarGlobal.showWindow)
                    arWaitbarGlobal.h = waitbar(j/n, arWaitbarGlobal.h, ...
                        sprintf('%i/%i   %2i%%   %s -> %s%s', ...
                        j, n, round(j/n*100), secToHMS(timeleft), ...
                        secToHMS(timeelapsed), funtext));
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



function hmstimestr = secToHMS(time)

days = floor(time/(60*60*24));
time = time - days*(60*60*24);

hours = floor(time/(60*60));
time = time - hours*(60*60);

minutes = floor(time/60);
time = time - minutes*60;

seconds = ceil(time);

% hmstimestr = sprintf('%id %ih %im %is', days, hours, minutes, seconds);

if(days>7)
    hmstimestr = sprintf('%id', days);
elseif(days>0)
    hmstimestr = sprintf('%id %ih', days, hours);
elseif(hours>0)
    hmstimestr = sprintf('%ih %im', hours, minutes);
elseif(minutes>0)
    hmstimestr = sprintf('%im %is', minutes, seconds);
else
    hmstimestr = sprintf('%is', seconds);
end
