arInit
arLoadModel('model_SEIR_InputSpline')
arLoadData('Germany__JHU_ge20_cum_fixS0')
arCompileAll

% This parameter is fixed and used to control at which time point the spline should be constant
arSetPars('tLast',max(ar.model.data.tExp),0,0,1,100,0); 
arAddEvent(1,1,max(ar.model.data.tExp));

arQplot('xy')
arLoadPars('fitted') % it's an illustration example, i.e. parameters were not checke for being necessarily the global optimum

arPlot

%%
x = linspace(-5,5,101);
y = 0.5*(1+(1-exp(-x))./(1+exp(-x)));
figure
plot(x,y,'LineWidth',2)
set(gca,'FontSize',14)
xlabel('spline input')
ylabel('b\_time\_dependence')
print -dpng SigmoidalFunction
