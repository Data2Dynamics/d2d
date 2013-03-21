
function hmstimestr = secToHMS(seconds)
hours = floor(seconds/3600);
seconds = seconds - hours*3600;
minutes = floor(seconds/60);
seconds = ceil(seconds - minutes*60);
hmstimestr = sprintf('%ih %im %is', hours, minutes, seconds);