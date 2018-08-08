% Illustration improvement of the retarded transient function to become
% time-Unit independent and to prevent numierical error raised by 10^t if
% t>300.
% 
% Example:
% t = linspace(0,10,101); 
% y = transFun2(t);
% y2 = transFun2(t*100,'timeUnitFactor',100);
% plot(t,y,t,y2)

function y = transFun2(t,varargin)
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

maxt = max(t);

timescale_sust = p.timescale_sust*p.timeUnitFactor/maxt;
timescale_trans = p.timescale_trans*p.timeUnitFactor/maxt;

tt = log10(10.^(10*t./maxt)+10.^(p.toffset_TF))-log10(1+10.^(p.toffset_TF));

y1 = p.amp_trans.*(1-exp(-tt./timescale_sust)).*exp(-tt./(timescale_sust+timescale_trans));
y2 = p.amp_sust.*(1-exp(-tt./timescale_sust));
y = y1+y2+p.offset_TF;
