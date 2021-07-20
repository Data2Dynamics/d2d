function [transf,ch_var,z_v] = new_var(V,ind,infi,rs_k,allVar,l_allVar,tmax,ode_n)
%%  FUNCTION NEW_VAR
%   Function that allows to make all the changes of variables. For it
% the function call "change.m" is necessary.
%
%   INPUT:      V --> Kernel
%               ind --> Number of columns of V
%               infi --> Infinitesimals
%               rs_k --> Coefficients of infinitesimal
%               allVar --> All variables
%               l_allVar --> allVar length
%               tmax --> Maximum degree of Lie series (passed to
%                        function "change.m")
%
%   OUTPUT:     transf --> Transformations vector
%               nuevas_variables --> Vector of the new variables
%%
    syms epsilon
    g_new_var=[];
    ch_var=[];
    for i=1:ind
        sv=V(:,i);
        aux1=subs(infi,rs_k,sv');
        a_g=gcd(aux1);
       if a_g~=1 && a_g~=0
           aux1=aux1./a_g;
       end
       aux1=simplify(aux1);
       g_new_var=[g_new_var;transpose(aux1)];
    end
    %   Remove combinations
    g_new_var=unique(g_new_var,'rows');
    s_g=size(g_new_var');
    R1=rref(g_new_var');
    r_R1=rank(R1);
    transf=[];
    R2=sort(R1,'descend');
    vec=zeros(s_g(1),1);
    vec(1)=1;
    %   transf vector
    for i=1:s_g(2)
       if R2(1:end,i)==vec
           transf=[transf;g_new_var(i,:)];
       end
    end
    %
    z_v=zeros(r_R1,l_allVar);
    for i=1:r_R1
        transf(i,:)=opt_vec(allVar,l_allVar,transf(i,:));
    end
    clear aux1
    for i=1:r_R1
        aux1=transf(i,:);
        if ode_n==0
            for j=1:l_allVar
                if aux1(j)~=0
                    [cc,tt]=coeffs(aux1(j),allVar(j));
                    %   See if it is a single term or are multiple
                    if length(cc)==1
                        [c3,c4]=numden(cc);
                        [cv,tv]=coeffs(c3,allVar);
                        if length(tv)==1 && abs(cv)==abs(cc) && abs(c4)==1
                            deg=feval(symengine,'degree',tt);              
                                if deg==0       %Translation
                                    nv=allVar(j)+cc*epsilon;
                                    z_v(i,j)=2;
                                elseif deg==1   %Scaling
                                    nv=exp(cc*epsilon)*allVar(j);
                                else            %Other
                                    nv=allVar(j)/(1-(deg-...
                                       1)*epsilon*allVar(j)^(deg-...
                                       1))^(1/(deg-1));
                                    z_v(i,j)=2;
                                end
                                aux1(j)=nv;                              
                        else
                            if exist('new')
                                aux1(j)=new(j);
                            else
                                new = change(allVar,l_allVar,...
                                    transf(i,:),tmax);
                                aux1(j)=new(j);
                            end 
                            z_v(i,j)=3;
                        end
                    else
                        if exist('new')
                            aux1(j)=new(j);
                        else
                            new = change(allVar,l_allVar,...
                                transf(i,:),tmax);
                            aux1(j)=new(j);
                        end
                        z_v(i,j)=3;
                    end
                else
                    aux1(j)=allVar(j);
                end
            end
            ch_var=[ch_var;aux1];
        else
           [new_var] =sol(allVar, aux1); 
           ch_var=[ch_var;transpose(new_var)];
        end

        clear new
    end
end

