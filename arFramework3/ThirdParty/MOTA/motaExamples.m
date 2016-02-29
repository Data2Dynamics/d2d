function K=motaExamples(Example)
% TEST - Examples
% motaExamples(Example) returns a (n X p)  matrix K as input for mota.m
% -------------------------------------------------------------------------

if Example==1
    dim=6;
    ll=200;
    K=5*(rand(ll,dim+1));
    K(:,1)=-K(:,2)+10+0.1*randn(ll,1);
    K(:,3)=5./(K(:,4).*K(:,5))+0.1*randn(ll,1);
    K(:,7)=0.1*ones(ll,1)+0.01*randn(ll,1);
end


if Example==2
    dim=5;
    ll=200;
    K=5*(rand(ll,dim+1));
    K(:,1)=3 * K(:,2) + 3 * K(:,3) + K(:,4) .* K(:,5)+2*randn(ll,1);
end


if Example==3
    dim=9;
    ll=200;
    K=5*(rand(ll,dim+1));
    K(:,1)=sqrt(1/50*(10*exp(K(:,2))+10*K(:,3).^3))+0.5*randn(ll,1);
    K(:,4)=5*K(:,5)+5*K(:,6)+randn(ll,1);
    K(:,7)=1./K(:,8)+0.3*K(:,9)+0.1*randn(ll,1);
end

if Example==4 
    dim=9;
    ll=200;
    K=5*(rand(ll,dim+1));
    K(:,1)=K(:,2).*K(:,3).*K(:,4)+0.5*randn(ll,1);
    K(:,5)=1/5*(3*K(:,6)+ 3*K(:,7)+K(:,8).*K(:,9))+0.1*randn(ll,1);
      
end

if Example==5
    dim=6;
    ll=200;
    K=5*(rand(ll,dim+1));
    K(:,1)=-K(:,2)+10+0.1*randn(ll,1);
    K(:,3)=5./(K(:,4).*K(:,5))+0.1*randn(ll,1);
    K(:,7)=0.1*ones(ll,1)+0.01*randn(ll,1);
end


if Example==6
    dim=9;
    ll=200;
    K=10*(rand(ll,dim+1));
    K(:,4)=15*(rand(ll,1));
    K(:,1) = K(:,2) + 0.2* K(:,3).*K(:,4)+0.5*randn(ll,1);
    K(:,6) = K(:,7)+eps*randn(ll,1);
    K(:,8) = 5*sin(K(:,9)) + 0.1*K(:,10).^2+0.5*randn(ll,1);
end


if Example==7
   dim=6;
   ll=200;
   K=5*rand(ll,dim+1);
   K(:,1)=-K(:,2) + 0.1*randn(ll,1);
   K(:,3)= 5./K(:,4) + 0.1*randn(ll,1);
   K(:,5)= 1/5*K(:,6).*K(:,7); 
end

if Example==8
   dim=3;
   ll=200;
   K=5*rand(ll,dim+1);
   K(:,1)=K(:,2).^2+ sin(K(:,3)) + 0.1*randn(200,1);  
   K(:,2)=K(:,2)+0.1*randn(200,1);
   K(:,3)=K(:,3)+0.1*randn(200,1);
   
end




