function [w, coh, phi, phi_low, phi_up] = kspec(x, y, dt, glbreite, wmax)
if (length(x) ~= length(y)) 
    error('x und y müssen gleiche länge haben!'); 
end
x = x(:);
y = y(:);
n = length(x);

%% Tapern
% tap = ones(1, n);
% tap = bartlett(n);
tap = hann(n);
tap = tap(:);
x = x.*tap;
y = y.*tap;

%% Fourieranalyse
if(exist('wmax', 'var'))
	[w, fxw] = fourier(x, dt, wmax);
	[w, fyw] = fourier(y, dt, wmax);
	else
	[w, fxw] = fourier(x, dt);
	[w, fyw] = fourier(y, dt);
end

%% Kreuzspektrum
csw = fxw.*conj(fyw);

%% Glätten (mit Hann Funktion(2*glbreite + 1))
if(~exist('glbreite', 'var'))
	glbreite = round(n/50);
end
csw = glaetten(csw, glbreite);
fxgl = glaetten(fxw.*conj(fxw), glbreite);
fygl = glaetten(fyw.*conj(fyw), glbreite);

coh = abs(csw) ./ sqrt(fxgl .* fygl);
phi = angle(csw);

%% Vertrauensintervall
%alpha = 0.05;
v = sum(hann(2*glbreite+1).^2) / 2;
tmp = sqrt((1./coh.^2 - 1) / v);
phi_low = phi - tmp;
phi_up = phi + tmp;