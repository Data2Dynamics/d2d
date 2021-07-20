function [A] = pol_sta_1(allVar,den_f,diffini,diff_den_f,diff_num_f,...
                        infi,l_f,l_in,num_f,rs_k)
%%  FUNCTION POL_OBS
%
%   Function for calculating the polynomial and system that must fullfied
%   ansatz constants for states to be to be fulfilled.
%   This function only applies to the 'uni' and 'partially' ansatz
%
%   INPUT:      allVar --> All variables
%               den_f --> f function denominator
%               diffini --> infi vector derivatives
%               diff_den_f --> den_f derivative
%               diff_num_f --> num_f derivative
%               infi --> Infinitesimals
%               l_f --> f function vector length
%               l_in --> infi vector length
%               num_f --> f function numerator
%               rs_k --> Coefficients of infinitesimal polynomial
%
%   OUTPUT:     A --> System matrix
%%
    for i=1:l_f
        nv=0;
        %   Polynomial
        for j=1:l_in
            nv=nv+infi(j)*(diff_num_f(i,j)*den_f(i)-...
                num_f(i)*diff_den_f(i,j));
        end
         eq_a=num_f(i)*den_f(i)*diffini(i)-nv; 

         %  System
        eq_c1=collect(eq_a,allVar);
        eq_ch1=children(eq_c1);
        l_eq_ch1=length(eq_ch1);
        eq_n1=[];
        for k=1:l_eq_ch1  
            [c_p1,~]=coeffs(eq_ch1(k),allVar);
            eq_n1=[eq_n1,c_p1];
        end
        
        %   Add to matrix
        [A1,~]=equationsToMatrix(eq_n1,rs_k);
        if i==1
            A=A1;
        else
            A=[A;A1];
        end
    end   
end

