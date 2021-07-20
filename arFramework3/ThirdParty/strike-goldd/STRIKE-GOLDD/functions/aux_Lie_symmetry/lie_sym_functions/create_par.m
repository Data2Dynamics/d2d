function [infi,rs_k] = create_par(allVar,m_aux,l_p,l_x,p,pMax,nw)
%%  FUNCION CREATE_PAR
% 
%   Function that creates partially varied infinitesimals
% 
%   INPUT:      allVar --> All variables
%               m_aux --> Auxiliar vector with the parameters and
%                         parametric ICS
%               l_p --> Parameter vector length
%               l_x --> State vector length
%               pMax --> Maximum infinitesimal degree
%               nw --> w vector length
%
%   OUTPUT:     infi --> Infinitesimal vector
%               rs_k --> Coefficients of infinitesimal polynomial
%
%%
    %   Creacion de la matriz con todas las combinaciones de parametros
    %   posibles (pueden superar el grado)
    rs_k=[];
    infi=[];
    for k=1:pMax-1    
        Pr1=m_aux.*(transpose(p));
        m_aux=[m_aux;nonzeros(Pr1)];
        m_aux=unique(m_aux);
    end
    vec_par=m_aux;
    %   Retirar del vector de combinaciones de parametros aquellas entradas
    %   que superar el maximo grado
    cont=1;
    while cont<=length(vec_par) 
       deg=feval(symengine,'degree',vec_par(cont));
       if deg>pMax
           vec_par(cont)=[];
       else
            cont=cont+1; 
       end
    end
    %   Creacion del vector de constantes r_s_d para los parametros
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

    %   Creacion de las potencias de los estados, uno a uno
    n_aux=[];
    for k=1:l_x+nw
        nv=[];
        for d=pMax+1:-1:1          
           s=d-1;
           nv=[nv,sym(allVar(k)^(s))];
        end
        nv=fliplr(nv);
        n_aux=[n_aux;nv];   
    end

    %   Multiplicacion del vector de las potencias de los estados
    %   por la matriz de parámetros
    s_m_aux=size(n_aux);
    infi2=[];
    l_vec_par=length(vec_par);
    for j=1:l_x+nw
        vec_sta=vec_par;

        %   Verificar que no supere el grado maximo y retirar aquellas 
        %   combinaciones que lo superen
        for i=2:s_m_aux(2)
              for k=1:l_vec_par
                deg=feval(symengine,'degree',vec_par(k));
                if deg+i-1<=pMax
                    vec_sta=[vec_sta;n_aux(j,i)*vec_par(k)];
                end
              end
        end
        vec_sta=unique(vec_sta);
        if isempty(vec_sta)==1
            vec_sta=n_aux; 
        end
        %   Creacion de las constantes r_s_d
        rs_a=[];
        for d=1:length(vec_sta)
            t=sym(['r_' num2str(j) '_'  num2str(d)]);
            rs_a=[rs_a,t];        
        end
        rs_k=[rs_k,rs_a];
        aux1=transpose(vec_sta)*transpose(rs_a);
        infi2=[infi2;aux1];   
    end
    infi=[infi2;infi];
end

