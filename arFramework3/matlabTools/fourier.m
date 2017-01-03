function [w, fw, wny] = fourier(data, dt, wmax)
n = length(data);

% Nyquistfrequenz
wny = pi/dt;
w = (0:n/2) * 2/n*wny;

if(exist('wmax','var'))
    if (wmax < wny) 
        w = w(w<wmax); 
    end
end

fw = fft(data);
fw = fw(1:length(w)) * 2 / sqrt(n);
