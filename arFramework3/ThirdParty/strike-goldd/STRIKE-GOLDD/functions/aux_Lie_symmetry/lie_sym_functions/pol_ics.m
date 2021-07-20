function [A] = pol_ics(allVar,den_ics,diff_den_ics,diff_num_ics,ics_n,...
                infi,known_ics,l_in,l_x,num_ics,rs_k,x)
%%  FUNCTION POL_ICS
%
%   Function for the creation of the polynomial and the system that must fullfied
%   ansatz coefficients for unknown ICS 
%
%   INPUT:      allVar --> All variables
%               den_ics --> ics_n denominator
%               diff_den_ics --> den_ics derivatives vector
%               diff_num_ics --> num_ics derivatives vector
%               ics_n --> ICS vector with infinitesimals coefficients
%               infi --> Infinitesimals
%               known_ics --> Known ICS
%               l_in --> infi length
%               l_x --> x length
%               num_ics --> ics_n numerator
%               rs_k --> Coefficients of infinitesimal polynomial
%               x --> States vector
%
%   OUTPUT:     A --> System matrix
%%
    ind_kics=find(known_ics);
    nv=subs(infi,x(ind_kics),ics_n(ind_kics));
    for i=1:l_x
        if known_ics(i)==1
            sum=0;
            %   Polynomial
            for j=1:l_in
                sum=sum+infi(j)*(diff_num_ics(i,j)*den_ics(i)-...
                   num_ics(i)*diff_den_ics(i,j)); 
            end
            sum=subs(sum,x,ics_n);

            %   System
            [eq_c_num,~]=numden(den_ics(i)^2*nv(i)-sum);
            eq_c_num=collect(eq_c_num,allVar);
            eq_ch=children(eq_c_num);
            l_eq_ch=length(eq_ch);
            eq_n=[];
            for k=1:l_eq_ch   
                [c_p,~]=coeffs(eq_ch(k),allVar);
                eq_n=[eq_n,c_p];
            end
            [A3,~]=equationsToMatrix(eq_n,rs_k);

            %   Add to matrix
            if exist('A','var')==0
                A=A3;
            else
                A=[A;A3];
            end
        end
    end
    if exist('A')==0
        A=[];
    end
end

