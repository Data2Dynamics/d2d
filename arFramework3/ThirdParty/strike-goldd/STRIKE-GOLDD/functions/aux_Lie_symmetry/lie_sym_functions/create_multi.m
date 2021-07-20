function [infi,rs_k] = create_multi(allVar,m_aux,l_p,l_x,p,pMax,nw,x)
%%  FUNCION CREATE_MULTI
% 
%   Function that creates multivariate infinitesimals.
% 
%   INPUT:      allVar --> All variables
%               m_aux --> Auxiliar vector with the parameters and
%                         parametric ICS
%               l_p --> Parameter vector length
%               l_x --> State vector length
%               pMax --> Maximum infinitesimal degree
%               nw --> w vector length 
%               x --> States vector
%
%   OUTPUT:     infi --> Infinitesimal vector
%               rs_k --> Coefficients of infinitesimal 
%
%%
    infi=[];
    rs_k=[];
    for k=1:pMax-1    
       Pr1=m_aux.*(transpose(p));
       m_aux=[m_aux;nonzeros(Pr1)];
       m_aux=unique(m_aux);
    end
    vec_par=m_aux;

    %   Remove from the vector all those combinations that exceed the
    %   maximum grade
    cont=1;
    while cont<=length(vec_par) 
       deg=feval(symengine,'degree',vec_par(cont));
       if deg>pMax
           vec_par(cont)=[];
       else
            cont=cont+1; 
       end
    end

    %   Creation of the r_s_d constants to be determined and the first
    % equations of the infinitesimals
    for k=1:l_p
       rs_a=[];
       for d=1:length(vec_par)
            t=sym(['r_' num2str(k+l_x+nw) '_'  num2str(d)]);
            rs_a=[rs_a,t];        
       end
       rs_k=[rs_k,rs_a];
       aux1=transpose(vec_par)*transpose(rs_a);
       infi=[infi;aux1];
    end

    %   Creation of the vector states from all combinations between them
    m_aux=[1;x];
    for k=1:pMax-1 
       Pr1=m_aux.*(transpose(x));   
       m_aux=[m_aux;nonzeros(Pr1)];
       m_aux=unique(m_aux);
    end
    s_m_aux=length(m_aux);
    l_vec_par=length(vec_par);
    for j=1:l_x
        vec_sta=vec_par;
        for i=1:s_m_aux
            deg_m=feval(symengine,'degree',m_aux(i));
              for k=1:l_vec_par
                deg=feval(symengine,'degree',vec_par(k));
                if deg+deg_m<=pMax
                    vec_sta=[vec_sta;m_aux(i)*vec_par(k)];
                end
              end
        end
        vec_sta=unique(vec_sta); 
    end
    if isempty(vec_sta)==1
        vec_sta=m_aux; 
    end 
    %   Unknown entries vector
    if nw~=0
        for k=l_x+1:l_x+nw
            nv=[];
            for d=pMax+1:-1:1          
               s=d-1;
               nv=[nv,sym(allVar(k)^(s))];
            end
            m_aux=[m_aux,nv];   
        end
        m_aux=unique(m_aux);
        s_m_aux=length(m_aux);
        l_vec_sta=length(vec_sta);
        for j=1+l_x:l_x+nw
            vec_w=vec_sta;
            for i=1:s_m_aux
                deg_m=feval(symengine,'degree',m_aux(i));
                  for k=1:l_vec_sta
                    deg=feval(symengine,'degree',vec_sta(k));
                    if deg+deg_m<=pMax
                        vec_w=[vec_w;m_aux(i)*vec_sta(k)];
                    end
                  end
            end
            vec_w=unique(vec_w); 
        end
        infi2=[];
        for i=l_x+1:l_x+nw
               rs_a=[];
               for d=1:length(vec_w)
                    t=sym(['r_' num2str(i) '_'  num2str(d)]);
                    rs_a=[rs_a,t];        
               end
               rs_k=[rs_k,rs_a];
               aux1=transpose(vec_w)*transpose(rs_a);
               infi2=[infi2;aux1];
        end
        infi=[infi2;infi];
    end

    %   Creation of r for states
    infi2=[];
    for i=1:l_x
           rs_a=[];
           for d=1:length(vec_sta)
                t=sym(['r_' num2str(i) '_'  num2str(d)]);
                rs_a=[rs_a,t];        
           end
           rs_k=[rs_k,rs_a];
           aux1=transpose(vec_sta)*transpose(rs_a);
           infi2=[infi2;aux1];
    end
    infi=[infi2;infi];
end

