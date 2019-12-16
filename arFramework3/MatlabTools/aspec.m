function [w, fw, sw, swgl, swgl_low, swgl_up] = aspec(data, dt, glbreite, wmax)
data = data(:);
n = length(data);

%% Tapern
% tap = ones(1, n);
% tap = bartlett(n);
tap = hann(n);
tap = tap(:);
data = data.*tap;

%% Fourieranalyse
if(exist('wmax','var'))
    [w, fw, wny] = fourier(data, dt, wmax);
else
    [w, fw, wny] = fourier(data, dt);
end

%% Autospektrum
sw = fw.*conj(fw);

%% Glätten (mit Hann Funktion(2*glbreite + 1))
if(~exist('glbreite','var'))
    glbreite = round(n/50);
end
swgl = glaetten(sw, glbreite);

%% Vertrauensintervall
alpha = 0.05;
v = sum(hann(2*glbreite+1).^2) / 2;
swgl_low = swgl / v * arChi2inv(alpha/2, v);
swgl_up = swgl / v * arChi2inv(1 - alpha/2, v);
