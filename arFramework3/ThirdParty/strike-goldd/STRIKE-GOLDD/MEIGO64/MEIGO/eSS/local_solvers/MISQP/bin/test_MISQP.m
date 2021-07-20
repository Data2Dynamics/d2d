                                      
xl=[0 0 0 0];                                          
xu=[10 10 10 10];                                    
x0=[3 4 5 1];                                               
nint=3;
ncont=1;
nbin=0;
m=3;
neq=3;
acc=1e-6;
FUN='test_f';
GRAD='evaluate_grad';

fobj='test_f';


run_misqp(ncont,nint,nbin,m,xl,xu,x0,acc,FUN,[],neq,fobj,[])

