function [A] = pol_obs(allVar,den_h,diff_den_h,diff_num_h,...
                        infi,l_h,l_in,num_h,rs_k)
%%  FUNCTION POL_OBS
%
%   Function for calculating the polynomial and system that must fullfied
%   ansatz constants for observations to be fulfilled
%
%   INPUT:      allVar --> All variables
%               den_h --> h function denominator
%               diff_den_h --> den_h derivative
%               diff_num_h --> num_h derivative
%               infi --> Infinitesimals
%               l_h --> h vector function length
%               l_in --> infi vector length
%               num_h --> h function numerator
%               rs_k --> Coefficients of infinitesimal polynomial
%
%   OUTPUT:     A --> System matrix
%%
    for i=1:l_h
        nv=0;
        %	Polynomial
        for j=1:l_in
           nv=nv+infi(j)*(diff_num_h(i,j)*den_h(i)-...
               num_h(i)*diff_den_h(i,j)); 
        end
        %   System
        eq_c=collect(nv,allVar);
        eq_ch=children(eq_c);
        l_eq_ch=length(eq_ch);
        eq_n=[];
        for k=1:l_eq_ch   
            [c_p,~]=coeffs(eq_ch(k),allVar);
            eq_n=[eq_n,c_p];
        end
        %   Add to matrix
        [A2,~]=equationsToMatrix(eq_n,rs_k);
        if i==1
            A=A2;
        else
            A=[A;A2];
        end
    end
end

