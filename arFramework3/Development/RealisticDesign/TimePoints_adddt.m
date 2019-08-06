function tT = TimePoints_adddt( p, changeT, folder_name)
%TIMEPOINTS Summary of this function goes here
%   p: Parameters on normal scale
global ar

%convertt = ar.model.data.convertt;
%dt = ar.model.data.dt;

if(~exist('p', 'var') || isempty(p))
    error('No parameters in TimePoint function')
end
%if(~exist('convertt', 'var') || isempty(convertt))
%    convertt = 1;
%end
if(~exist('changeT', 'var') || isempty(changeT))
    changeT = 1;
end
%if(~exist('dt', 'var') || isempty(dt))
%    error('No time shift dt given. Use function TimePoints instead of TimePoints_adddt, or pass time shift dt in function declaration.')
%end
if(~exist('folder_name', 'var') || isempty(folder_name))
    folder_name = '';
end

if any(any(isnan(p)))==1
error 'A Fitting did not work, a parameter is NaN'
end
if any(any(isinf(p)))==1
v1 = find(isinf(p));
error('Observable %i is constant. No time point estimation reasonable.',v1)
end

% norm
T = ( 0.2 + p(:,2)./(p(:,1)+p(:,2)) ) * (p(:,1)+p(:,2));
n = ceil(10.^( (0.7*p(:,2)./(p(:,1)+p(:,2))) * (p(:,1)+p(:,2)) ));

% halflog
%T = 37.4 + 21.9*p(:,1) + 142.4*p(:,2);
%n = ceil(12.9 - 0.4*p(:,1) - 0.6*p(:,2));
%add
%     T = 1.85 + 0.15 *(p(:,1)+p(:,2));
%     n  = ceil(10.^(1.12 -0.03*(p(:,1)+p(:,2))));
%sep
%     T = 1.80 + 0.11 *p(:,1) + 0.21 *p(:,2);
%     n  = ceil(10.^(   1.06 -0.09 *p(:,1) +0.05 *p(:,2)   ));
% multi
%    T = 1.76 + 0.15*p(:,1) + 0.24*p(:,2) -0.03*p(:,1).*p(:,2);
%    n = ceil(10.^(1.04 -0.07*p(:,1) + 0.06*p(:,2) -0.01*p(:,1).*p(:,2)));    
%    T = log10(10.^(T)*changeT);

tT=nan(max(n),length(T));
for k=1:length(T)
    tT(1:ceil(n(k)),k)=logspace(0,T(k),n(k))-1;
    %tT(:,k) = tT(:,k)./convertt(k)+dt(k);
end

tT = round(tT,2)
tT(tT>3) = round(tT(tT>3)/5,1)*5;
tT(tT>30) = round(tT(tT>30)/50,1)*50;
tT(tT>300) = round(tT(tT>300)/500,1)*500;
tT(tT>3000) = round(tT(tT>3000)/5000,1)*5000;

if changeT  == 1
    changeT = [];
end

ar.model.data.tExp = tT;
xlswrite(['RealisticDesign' num2str(round(changeT)) '/TimePoints.xls'],tT);
fprintf('Realistic time Points are assigned.\n');    
    
end

