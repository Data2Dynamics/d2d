function [A] = pol_sta_2(allVar,den_f,diffini,diff_den_f,diff_num_f,...
                        infi,l_f,l_in,num_f,rs_k)
%%  FUNCTION POL_STA_2
%
%   Function for calculating the polynomial and system that must fullfied
%   ansatz constants for states to be met.
%   This function only applies to the 'multi' ansatz
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
    for k=1:l_f
        % 1st polynomial || Construction
        nv=0;
        nm=0;
        for i=1:l_in
            nv=nv+infi(i)*(diff_num_f(k,i)*den_f(k)-...
                num_f(k)*diff_den_f(k,i))*prod(den_f(1:end ~=k));
        end
        for j=1:l_f
            nm=nm+diffini(k,j)*den_f(k)*num_f(j)*prod(den_f(1:end ~=j));      
        end
         eq_a=nm-nv;
         
        % System
        eq_c1=collect(eq_a,allVar);
        eq_ch1=children(eq_c1);
        l_eq_ch1=length(eq_ch1);
        eq_n1=[];
        for l=1:l_eq_ch1  
            [c_p1,~]=coeffs(eq_ch1(l),allVar);
            eq_n1=[eq_n1,c_p1];
        end
        
        % Add to matrix
        [A1,~]=equationsToMatrix(eq_n1,rs_k);
        if k==1
            A=A1;
        else
            A=[A;A1];
        end
    end 
end

