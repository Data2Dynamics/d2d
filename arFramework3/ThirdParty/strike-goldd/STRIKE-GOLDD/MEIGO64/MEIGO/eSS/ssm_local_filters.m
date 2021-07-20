% BUG FIXED JUNE-2014 - JRB

if stage_1==1
    %                         Aplicamos el solver local sobre el mejor
    %                         de los hijos, incluidos los del Refset.
    [minimo,I]=min(child_values_penalty);
    Threshold=min(minimo,fbest);
    if minimo<fbest
        x0=child(I,:);
        f0=minimo;
    else
        x0=xbest;
        f0=fbest;
    end
    apto=2;
    %Pasamos al stage 2
    n_critico=local_n2;
    %Llamamos al solver local

elseif stage_2==1;
    if local_bestx
        apto=2;
        if use_local
            x0=xbest;
            f0=fbest;
        else
            apto=0;
        end
    else
        %Tenemos que ver si pasa los filtros
        %De momento decimos que el punto no pasa
        %los filtros
        apto=0;

        %Si el filtro de merito esta activado
        if merit_filter
            if child_values_penalty(k)<Threshold
                wait_Threshold=0;
                apto=apto+1;
            else
                wait_Threshold=wait_Threshold+1;
                if wait_Threshold==wait_th_limit;
                    Threshold=Threshold+thfactor*(1.0+abs(Threshold));
                    wait_Threshold=0;
                end
            end
            %Si no esta activado, directamente pasa
        else
            apto=apto+1;
        end

        %Asi hemos hecho independientes el filtro de
        %merito y el de distancia. Antes los tenia como
        %dependientes. Si no pasaba el de merito no
        %evaluaba el de distancia


        %Vamos a por el filtro de la distancia
        if distance_filter
            lejano=0;
            for m=1:size(local_solutions,1)
                %Hallamos la distancia entre el punto y las soluciones
                separacion=norm((child(k,:)-local_solutions(m,1))./(x_U-x_L));
                %Si la separacion al minimo es mayor que maxdist para
                %ese minimo
                if  separacion>maxdist(m)
                    lejano=lejano+1;
                    %Reseteamos la espera para ese
                    %minimo
                    wait_maxdist(m)=0;
                else
                    wait_maxdist(m)=wait_maxdist(m)+1;
                    %disp('NO pasa filtro de la distancia')
                    if wait_maxdist(m)==wait_maxdist_limit;
                        maxdist(m)=maxdist(m)*(1-maxdistfactor);
                        wait_maxdist(m)=0;
                    end
                    break
                end
            end
            if lejano==size(local_solutions,1)
                %disp('Pasa filtro de la distancia')
                apto=apto+1;
            end
            %Si no esta activado, pasamos de el
        else
            apto=apto+1;
        end
        %Vemos si esa solucion la tenemos
        %ya porque podria ser que
        %iniciaramos el local solver con x0
        %y ya tuvieramos alguna
        if apto==2
            x0=child(k,:);
            f0=child_values_penalty(k);
        end
    end
end
if apto==2
    n_minimo=0;

    if iterprint
        fprintf('Call local solver: %s \n', upper(local_solver))
        fprintf('Initial point function value: %f \n',f0);

        tic
    end
    if isFjacDefined
        [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,fjac,...
            local_solver,local_iterprint,local_tol,weight,nconst,tolc,opts.local,varargin{:});
    else
        [x,fval,exitflag,numeval]=ssm_localsolver(x0,x_L,x_U,c_L,c_U,neq,ndata,int_var,bin_var,fobj,[],...
            local_solver,local_iterprint,local_tol,weight,nconst,tolc,opts.local,varargin{:});
    end
    if iterprint
        fprintf('Local solution function value: %f \n',fval);
        fprintf('Number of function evaluations in the local search: %i \n',numeval);
        fprintf('CPU Time of the local search: %f  seconds \n\n',toc);
    end

    %Ojo, que solnp no da el numero de evaluaciones
    %de funcion

    nfuneval=nfuneval+numeval;

% BUG FIXED JUNE-2014 - JRB - IF BLOCK BELOW REMOVED
%    if stage_1==1
%        stage_1=0;
%        stage_2=1;
%    end

    %Evaluate the local solution
    [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(x,x_L,x_U,fobj,nconst,c_L,c_U,...
        tolc,weight,int_var,bin_var,varargin{:});
    nfuneval=nfuneval+1;
    if include & pena<=tolc
	% BUG FIXED JUNE-2014 - JRB INCLUDING IF BLOCK LINE BELOW
	if stage_1==1, stage_1=0; stage_2=1; end
        if val_penalty<fbest
            fbest=val_penalty;
            xbest=x;
            if fbest<fbest_lastiter
                Results.f=[Results.f fbest];
                Results.x=[Results.x; xbest];
                Results.time=[Results.time cputime-cpu_time];
                Results.neval=[Results.neval nfuneval];
                if plot_results==1
                    stairs(Results.time,Results.f)
                    xlabel('CPU Time (s)');
                    ylabel('Objective Function Value (-)');
                    title(strcat('Best f:  ',num2str(fbest)))
                    drawnow
                end
            end
        end



        %Chequeamos el criterio de parada
        if nfuneval>= maxeval
            fin=1;
            blacklist=ones(dim_refset);
            break
        elseif cputime-cpu_time>=maxtime
            fin=2;
            blacklist=ones(dim_refset);
            break
        elseif not(isempty(vtr))
            if fbest<=vtr
                fin=3;
                blacklist=ones(dim_refset);
                break
            end
        end





        %Añadimos ese valor a los hijos

        child=[child;x];

        if nconst
            [fffff,ggggg]=feval(fobj,x,varargin{:});
        else
            fffff=feval(fobj,x,varargin{:});
            ggggg=0;
        end

        child_values=[child_values;fffff];

        if size(ggggg,2)==1
            ggggg=ggggg';
        end

        child_nlc=[child_nlc;ggggg];

        nfuneval=nfuneval+1;

        child_penalty=[child_penalty;ssm_penalty_function(x,child_nlc(k,:),c_L,c_U,tolc)];
        child_values_penalty=[child_values_penalty;child_values(k)+weight*child_penalty(k)];


        %Añadimos ese valor a los hijos. Habra que
        %añadir la opcion de que esto no sea asi
        %Sustituimos a su padre

        %En posteriores versiones habra que ver
        %si compensa que el hijo mejorado
        %reemplace al padre o si por el
        %contrario se añade para que se
        %combinen juntos.

        %El nuevo Threshold pasa a ser el minimo
        %obtenido.
        if fval<Threshold
            Threshold=fval;
        end

        if not(isempty(local_solutions))
            %Vamos a ver si esa solucion ya la
            %teniamos
            %En principio decimos que es un nuevo
            %local
            adicionar_local=1;
            %Comprobamos si ese minimo ya lo
            %tenemos
            [f,ind]=ssm_isdif2(x,local_solutions,1e-2,1);
            %Si esta cerca de alguno
            %actualizamos las distancias
            if f
                for i=length(ind)
                    %Si es el mismo minimo
                    %Actualizamos el maxdist para ese minimo
                    maxdist(ind(i))=norm((child(k,:)-local_solutions(ind(i),:))./(x_U-x_L));
                    %Y no añadimos a la lista de
                    %minimos puesto que ya lo
                    %tenemos
                    adicionar_local=0;
                    %Una vez sabemos que es uno de
                    %los que tenemos, ya podemos
                    %salir del bucle for. No va a
                    %ser ninguno de los otros si el
                    %filtro esta bien hecho
                    %break
                end
            end

        else
            %Si no habia ninguna solucion local, se
            %añade
            adicionar_local=1;
        end
        %Si tenemos que añadir a la lista de
        %soluciones locales, lo hacemos
        if adicionar_local==1
            local_solutions=[local_solutions;x];
            local_solutions_values=[local_solutions_values;fval];
            maxdist=[maxdist;norm((child(k,:)-x)./(x_U-x_L))];
            wait_maxdist=[wait_maxdist; 0];
        end
        %Tanto si hemos añadido minimos locales
        %como si hemos actualizado las maxdist
        %Vamos a chequear si las maxdist se
        %solapan, siempre y cuando el filtro de la
        %distancia este activado, porque si no no
        %tiene sentido
        if distance_filter
            nlocals=size(local_solutions,1);
            if nlocals>1
                ncomb_local=nchoosek(1:nlocals,2);
                number_of_comb=(nlocals^2-nlocals)/2;
                for i=1:number_of_comb
                    index=ncomb_local(i,:);
                    ecart=norm((local_solutions(index(1),:)-local_solutions(index(2),:))./(x_U-x_L));
                    while ecart<maxdist(index(1))+maxdist(index(2))
                        maxdist(index(1))=0.9*maxdist(index(1));
                        maxdist(index(2))=0.9*maxdist(index(2));
                    end
                end
            end
        end
    end
end