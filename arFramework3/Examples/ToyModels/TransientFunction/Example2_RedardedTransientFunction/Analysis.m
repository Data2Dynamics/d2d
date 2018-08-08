%% Illustration of the transient function with toffset (for modelling retarded transient response)
t = linspace(0,10,101);
t_offsets = [-1,0,1,2,5];

h = [];
legstr = cell(0);
clf
hold on
for i=1:length(t_offsets)
    toffset_TF = t_offsets(i);
    h(end+1)=plot(t,log10(10.^t+10^toffset_TF)-log10(1+10^toffset_TF));
    legstr{end+1} = ['toffset\_TF=',num2str(toffset_TF)];    
end

for i=1:length(t_offsets)
    toffset_TF = t_offsets(i);
    plot(t,transFun(t,toffset_TF),'Color',get(h(i),'Color'),'LineWidth',2);
end
title(strrep('log_{10}(10^t+10^{toffset_TF})-log_{10}(1+10^{toffset_TF})','_','\_'))
xlabel('t')
legend(legstr{:},'Location','NorthWest')
print -dpng TimeOffset_TransientFunction


%% This plot illustrates 
% a) the numerical problem if times >300 occur
% b) the problem, that the smoothness of the time_transformation is not
% unit-independent
t = linspace(0,10,101); 
y = transFun(t);
y2 = transFun(t*100,'timeUnitFactor',100);
plot(t,y,t,y2)

title(strrep('log_{10}(10^t+10^{toffset_TF})-log_{10}(1+10^{toffset_TF})','_','\_'))
xlabel('t')
legend('tUnitFactor=1','tUnitFactor=100','Location','NorthEast')
print -dpng TimeOffset_TransientFunction_timeUnitProblem

%% This time_transformation 
% unit-independent
t = linspace(0,10,101); 
y = transFun2(t);
y2 = transFun2(t*100,'timeUnitFactor',100);
plot(t,y,'-',t,y2,'--')

title(strrep('log_{10}(10^t+10^{toffset_TF})-log_{10}(1+10^{toffset_TF})','_','\_'))
xlabel('t')
legend('tUnitFactor=1','tUnitFactor=100','Location','NorthEast')
print -dpng TimeOffset_TransientFunction_timeUnitInvariant

%%
t = linspace(0,10,101); 
t_offsets = [-1,0,1,2,5];
tfacs = logspace(0,3,4);
close all
for j=1:length(tfacs)
    tfac = tfacs(j);
    legstr = cell(0);
    subplot(2,2,j)
    hold on
    for i=1:length(t_offsets)
        toffset_TF = t_offsets(i);
        plot(t,transFun2(t*tfac,'timeUnitFactor',tfac,'toffset_TF',toffset_TF),'LineWidth',2);
        legstr{end+1} = ['t\_offset=',num2str(toffset_TF)];
    end
    legend(legstr{:});
    title(['timeUnitFactor = ',num2str(tfac)])
    xlabel('t [timeUnitFactor]');
end


%% Solution 1: Taylor expansion of 10^t
clear t
syms t
T = taylor(10^t,t)

t = linspace(0,10,101);
toffset_TF = 1;
close all
plot(t,log10(10.^t+10^toffset_TF)-log10(1+10^toffset_TF))
hold on
for order=1:5
    plot(t,log10(TaylorApprox(t,order)+TaylorApprox(toffset_TF,order))-log10(1+TaylorApprox(toffset_TF,order)));
end
legend('Exact','order=1','order=2','order=3','order=4','order=5');
title('Does not work!');


%% Solution 1b: Taylor expansion of the whole transformation
clear t
syms t toff
T = taylor(log10(10^t+10^toff)-log10(1+10^toff))




