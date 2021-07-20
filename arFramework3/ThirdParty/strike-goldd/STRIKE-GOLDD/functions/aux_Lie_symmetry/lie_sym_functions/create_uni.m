function [infi,rs_k] = create_uni(allVar,l_allVar,pMax)
%%  FUNCION CREATE_UNI
% 
%   Function that creates univariate infinitesimals
% 
%   INPUT:      allVar --> All variables
%               l_allVar --> allVar length
%               pMax --> Maximum infinitesimal degree
%
%   OUTPUT:     infi --> Infinitesimal vector
%               rs_k --> Coefficients of infinitesimal polynomial 
%
%%
    infi=[];
    rs_k=[];   
    for k=1:l_allVar
       temp=0;
           for d=1:pMax+1
              s=d-1;
              t=sym(['r_' num2str(k) '_'  num2str(s)]);
              rs_k=[rs_k,t]; 
              temp=temp+sym(rs_k(end)*sym(allVar(k)^s));
           end
       infi=[infi;temp];
    end
end

