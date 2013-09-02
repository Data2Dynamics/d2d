
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