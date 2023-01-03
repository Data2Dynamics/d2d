function tT = arSetTimePoints( Table, nx, ny, convertt, folder_name)

global ar

if(~exist('Table', 'var') || isempty(Table))
    error('No parameters in TimePoint function')
end
if(~exist('convertt', 'var') || isempty(convertt))
    convertt = 1;
end
if(~exist('folder_name', 'var') || isempty(folder_name))
    folder_name = '';
end

sust_times_trans = Table.timescale_sust.*Table.timescale_trans;
max_sus_model = nanmax(Table.timescale_sust);
max_tr_model = nanmax(Table.timescale_trans);
min_sus_model = nanmin(Table.timescale_sust);
min_tr_model = nanmin(Table.timescale_trans);

% aus 2017
%T = 1.76 + 0.15*tsus + 0.24*ttrans -0.03*tsus.*ttrans+0.01*toff;
%n = ceil(10.^(1.04 -0.07*tsus + 0.06*ttrans -0.01*tsus.*ttrans));
% 2022:
T =      10.^(2.17 +0.00034*Table.toffset_TF +0.0035/nx +0.014/ny +0.011*Table.amp_sust -0.022*Table.amp_trans +0.67*min_sus_model +0.08*min_tr_model +0.03*Table.timescale_sust +0.143*Table.timescale_trans +0.03*sust_times_trans +randn(size(Table,1),1)*0.03); % 0.3
n = round(10.^(1.04 -0.00032*Table.toffset_TF +0.0015/nx -0.016/ny +0.013*Table.amp_trans -0.07*max_sus_model +0.07*max_tr_model +0.049*Table.timescale_sust +0.05*Table.timescale_trans -0.033*sust_times_trans +randn(size(Table,1),1)*0.02)); % 0.2
lambda =  0.667 -0.00017*Table.toffset_TF -0.0038/nx +0.0039/ny -0.013*Table.amp_trans -0.121*max_sus_model +0.086*max_tr_model -0.04*Table.timescale_trans      +randn(size(Table,1),1)*0.01;
potenz =  1.75 +0.00022*Table.toffset_TF +0.0086/nx -0.009/ny +0.036*Table.amp_trans +0.14*max_sus_model -0.18*max_tr_model +0.08*Table.timescale_sust +0.14*Table.timescale_trans -0.032*sust_times_trans +randn(size(Table,1),1)*0.03;
% nur observable variable, not model parameter
T =      10.^(1.58 +0.00151*Table.toffset_TF -0.0046/nx +0.016/ny -0.019*Table.amp_trans +0.13*Table.timescale_sust +0.28*Table.timescale_trans +0.042*sust_times_trans +randn(size(Table,1),1)*0.04); % 0.3
n = round(10.^(1.02 -0.00029*Table.toffset_TF +0.0016/nx -0.014/ny +0.015*Table.amp_trans +0.037*Table.timescale_sust +0.059*Table.timescale_trans -0.03*sust_times_trans +randn(size(Table,1),1)*0.02)); % 0.2
lambda =  0.616 -0.00015*Table.toffset_TF -0.0039/nx +0.0063/ny -0.006*Table.amp_sust -0.006*Table.amp_trans -0.014*Table.timescale_sust -0.044*Table.timescale_trans     +randn(size(Table,1),1)*0.01;
potenz =  1.73 +0.0073/nx -0.013/ny +0.031*Table.amp_trans +0.08*Table.timescale_sust +0.098*Table.timescale_trans -0.035*sust_times_trans +randn(size(Table,1),1)*0.04;

tT=nan(max(n),length(T));
for k = 1:length(T)
    x = linspace(0,1,n(k));
    tT(1:n(k),k) = (lambda(k).*x+(1-lambda(k)).*x.*((exp(log(2)*x)-1).^potenz(k))).*T(k);
end

tT = tT./convertt'; % biomodel was scaled to range(t)>10 || <100

%tT = round(tT,2,'significant');
mag = 10.^(8:-1:-8);
for i=1:length(mag)
    tT(tT>mag(i)) = round(tT(tT>mag(i))*2,-log10(mag(i)))/2;
end
tT

ar.model.data.tExp = tT;
writematrix(tT,'RealisticDesign/TimePoints.txt');
%writematrix(tT2,'RealisticDesign/TimePoints2.txt','WriteMode','overwrite');
fprintf('Realistic time Points are assigned.\n');    
    
end

