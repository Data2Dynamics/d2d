function [new_vec] = opt_vec(allVar,l_allVar,transf)
%%  FUNCTION OPT_VEC
%
%   Function that looks for the vector of transformations with less changes in
%   parameters that depend on Lie series.
%
%   INPUT:      allVar --> All variables
%               l_allVar --> allVar length
%               transf --> Infinitesimal generator vector
%
%   OUTPUT:     
%               new_vec --> Vector with the new infinitesimal generators
%
%%
    z_v=zeros(1,l_allVar);
    for j=1:l_allVar
        if transf(j)~=0
            [cc,~]=coeffs(transf(j),allVar(j));
            if length(cc)==1
                vec2=transf./cc;
                i=find(abs(cc)==abs(allVar));
                for l=1:l_allVar
                    if vec2(l)~=0
                        [c3,c4]=numden(vec2(l));
                        [~,tv]=coeffs(c3,allVar);
                        if length(tv)==1 && tv==allVar(l) && abs(c4)==1
                            z_v(i)=z_v(i)+1;
                        end
                    end
                end
            end
        end
    end
    [m,i]=max(z_v);
    if m~=0
        vec_f=transf./allVar(i);
    else
        vec_f=transf;
    end
    new_vec=vec_f;
end

