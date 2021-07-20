%             [NAGWare Gateway Generator]
%
%Copyright (c) 1993-97 by the Numerical Algorithms Group Ltd 2.0a
%
%n2fb
%
%n                                     integer
%p                                     integer
%x (p)                                 real
%b (2,p)                               real
%iv (liv)                              integer
%liv                                   integer
%lv                                    integer
%v (lv)                                real
%ui (:)                                integer
%ur (:)                                real
%
%[x,iv,v,ur] = n2fb(n,p,x,b,iv,liv,lv,v,ui,ur)
%
%
% MRF 14/06/2005

function fobj = call_n2fb(X,N,P,xl,xu,local_solver)


%tic;

% *** THE FOLLOWING CODE IS FOR CALLING N2FB.
%
  %   N = 21*8*16;
  %   P = 36;
     LIV = 82+4*P;
     LTY = N;
     LV = 98+P*(3*P+25)/2+7+P+N*(P+2)+P*(P+15)/2;
%
% *** Y VALUES...
%
    for i = 1:N
        TY(i,2) = 0;
    end
%
% ***  T VALUES...
%
    for i = 1:N
        TY(i,1) = 0;
    end
%
% *** SUPPLY LEAD DIMENSION OF TY IN UI(1)...
% *** (MOST COMPILERS WOULD LET US SIMPLY PASS LTY FOR UI,
% *** BUT SOME, E.G. WATFIV, WILL NOT.)
%
      UI(1) = LTY;
      IV(1:LIV) = 0;
      IV(17) = 2000; % MXFCAL
      IV(18) = 3000; % MXITER
      V(1:LV) = 0;


      % *** SUPPLY BOUNDS...
B(1,:)=xl;
B(2,:)=xu;



%Call the MEX function
%
if strmatch(local_solver,'n2fb')
    [X,IV,V,TY] = wn2fb(N,P,X,B,IV,LIV,LV,V,UI,TY);
else
    %dn2fb
    [X,IV,V,TY] = wdn2fb(N,P,X,B,IV,LIV,LV,V,UI,TY);
end

fobj=X';
% CPU_time=toc
