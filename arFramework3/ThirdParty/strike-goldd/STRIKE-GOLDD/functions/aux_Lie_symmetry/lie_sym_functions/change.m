function [nv] = change(allVar,l_allVar,transf,tmax)
%%  FUNCTION CHANGE
%   It finds the change of variable when the transformation
%   depends on two or more parameters or, in case the generator
%   infinitesimal does not correspond to the variable.
%   INPUT:      allVar --> All variables
%               transf --> Transformations vector
%               tmax   --> Maximum degree for Lie series
%
%   OUTPUT:     new --> New variables
%%
    syms epsilon
    k_ind=find(transf);
    transf=nonzeros(transf);
    vec_it=transf;
    transf=transpose(transf);
    n_allVar=allVar(k_ind);
    l_nallVar=length(n_allVar);
%   Initialization for grade 0 and 1
    new=n_allVar;
    new=new+epsilon.*vec_it;
%   Calculation of higher degrees
    for n=2:tmax
       %    Powers of the generator vector 
        A_it=[];
        for i=1:l_nallVar
           aux=[];
           for j=1:l_nallVar
               aux=[aux,diff(vec_it(j),n_allVar(i))];
           end
           A_it=[A_it;aux];
        end
        vec_it=(transf)*A_it;
        vec_it=transpose(vec_it);
        vec_it=simplify(vec_it);
        %   Product of the generator vector and variables 
        new=new+((epsilon^n/factorial(n)).*vec_it);
    end
    nv=allVar;
    nv(k_ind)=new;
end

