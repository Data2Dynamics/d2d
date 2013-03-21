function out = glaetten(in, glbreite)

xsize = length(in);
out = zeros(size(in));

% %% Mit Rechteck glaetten (schlecht)
% for j=1:xsize
% 	glstart = j-glbreite;
% 	glstop = j+glbreite;
% 	if(glstart < 1) glstart = 1; end
% 	if(glstop > xsize) glstop = xsize; end
% 	out(j) = sum(in(glstart:glstop)) / (glstop - glstart + 1);
% end

%% Mit Hann Funktion glaetten (gut)
for j=1:xsize
    glstart = j-glbreite;
    glstop = j+glbreite;
    window = my_hann(2*glbreite + 1);
    if(glstart < 1)
        window = window(abs(glstart)+2:2*glbreite+1);
        glstart = 1;
    end
    if(glstop > xsize)
        window = window(1:2*glbreite+1-(glstop-xsize));
        glstop = xsize;
    end
    tmpin = in(glstart:glstop);
    tmpin = tmpin(:);
    window = window(:);
    out(j) = sum(tmpin.*window) / sum(window);
end

function y=my_hann(xsize)
if(xsize==2)
    y=[0.5 0.5];
elseif(xsize==1)
    y=1;
elseif(xsize<1)
    error('my_hann(size>0) !!!');
else
    x = -(xsize-1)/2:(xsize-1)/2;
    y = 0.5*(1+cos(2*pi*x/(xsize-1)));
end


