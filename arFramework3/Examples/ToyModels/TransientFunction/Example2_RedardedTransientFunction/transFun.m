% % Illustration of the retarded transient function, i.e. the classical
% transient function with "time_offset"
% 
% Example:
% toffset_TF = 1; 
% plot(t,transFun(t))

function y = transFun(t,varargin)
p  = inputParser;
addRequired(p,'t',@isnumeric);
addOptional(p,'toffset_TF',0,@isnumeric);
addOptional(p,'timeUnitFactor',1,@isnumeric);
addOptional(p,'amp_sust',1,@isnumeric);
addOptional(p,'amp_trans',4,@isnumeric);
addOptional(p,'timescale_sust',1,@isnumeric);
addOptional(p,'timescale_trans',1,@isnumeric);
addOptional(p,'offset_TF',0,@isnumeric);
parse(p,t,varargin{:});

p = struct(p.Results);

timescale_sust = p.timescale_sust*p.timeUnitFactor;
timescale_trans = p.timescale_trans*p.timeUnitFactor;

tt = log10(10.^t+10.^p.toffset_TF)-log10(1+10.^p.toffset_TF);

y1 = p.amp_trans.*(1-exp(-tt./timescale_sust)).*exp(-tt./(timescale_sust+timescale_trans));
y2 = p.amp_sust.*(1-exp(-tt./timescale_sust));
y = y1+y2+p.offset_TF;
