
InputSpline = linspace(-5,5,101);

y1 = 0.5*(1+(1-exp(-InputSpline))./(1+exp(-InputSpline)));

y2 = 1./(1+exp(-InputSpline));
y3 = exp(InputSpline)./(1+exp(InputSpline));

plot(InputSpline,y1,'-',InputSpline,y2,'--',InputSpline,y3,':')
legend('Parametrization 1','Parametrization 2','Parametrization 3')


