function [xout,fout,improve,evals_local, Results]=rvnds_local(x0,f0,fobj,x_L,x_U,search_type,shaked_vars,decomp,nvar,maxevals,maxtime,cpu_time,fbest,xbest,nfuneval, Results, varargin)

evals_local=0;
last_var_used=0;

fin=0;

xcurrloc=x0;
fcurrloc=f0;

if decomp
    u_vars=shaked_vars;     
    n_u_vars=numel(shaked_vars);
else
    u_vars=1:nvar;
    n_u_vars=nvar;
end

improve_solution=0;

while ~fin
    x_improvements=[];
    f_improvements=[];
    vars_improving=[];
    exit_local=0;

    possible_directions=repmat([-1 1],n_u_vars,1);


    %Find variables touching bounds
    aaa=find(xcurrloc(u_vars)==x_L(u_vars));
    bbb=find(xcurrloc(u_vars)==x_U(u_vars));

    %Adjust the possible directions. i.e., variables touching bounds
    %can only go in one direction
    possible_directions(aaa,1)=0;
    possible_directions(bbb,2)=0;

    %Choose the order of the variables to explore randomly
    ccc=randperm(n_u_vars);
    possible_directions=possible_directions(ccc,:);

    local_improvement=0;
    for i=1:n_u_vars
        if u_vars(ccc(i))~=last_var_used
            for j=1:2   %2 possible directions: +1 and -1 in the best case
                exit_local=0;       %This has to be here or 1 level up???
                if possible_directions(i,j)     %If we can move to that direction
                    xtemp=xcurrloc;
                    xtemp(u_vars(ccc(i)))=xcurrloc(u_vars(ccc(i)))+possible_directions(i,j);
                    ftemp=feval(fobj,xtemp, varargin{:});
                    evals_local=evals_local+1;
					
                    if ftemp<fcurrloc
						local_improvement=1;
						improve_solution=max(improve_solution,1);  %Avoid improve_solution=1 if it was =2 previously. This can only happen in best improvement

                        if ftemp<fbest;
                            fbest=ftemp;
                            xbest=xtemp;
                            improve_solution=2;
							
							Results.xbest=xbest;
							Results.fbest=fbest;
							Results.func=[Results.func; fbest];
							Results.x=[Results.x;xbest];
							Results.time=[Results.time; cputime-cpu_time];
							Results.neval=[Results.neval; nfuneval+evals_local];
							Results.numeval=nfuneval+evals_local;
							Results.cpu_time=cputime-cpu_time;
							
							fprintf('NFunEvals: %i  Bestf: %g      CPUTime: %f   \n',nfuneval+evals_local,fbest,cputime-cpu_time);
                        end



                        %%go beyond
                        parar=0;
                        while ~parar
                            if xtemp(u_vars(ccc(i)))==x_L(u_vars(ccc(i))) | xtemp(u_vars(ccc(i)))==x_U(u_vars(ccc(i)))
                                %touching bounds!!!
                                parar=1;
                            else
                                %Continue searching in the same direction
                                xtemp2=xtemp;
                                xtemp2(u_vars(ccc(i)))=xtemp(u_vars(ccc(i)))+possible_directions(i,j);
                                ftemp2=feval(fobj,xtemp2,varargin{:});
                                evals_local=evals_local+1;
                                if ftemp2<fbest;
                                    fbest=ftemp2;
                                    xbest=xtemp2;
                                    improve_solution=2;
									Results.xbest=xbest;
									Results.fbest=fbest;
									Results.func=[Results.func; fbest];
									Results.x=[Results.x;xbest];
									Results.time=[Results.time; cputime-cpu_time];
									Results.neval=[Results.neval; nfuneval+evals_local];
									Results.numeval=nfuneval+evals_local;
									Results.cpu_time=cputime-cpu_time;
                                    fprintf('NFunEvals: %i  Bestf: %g      CPUTime: %f   \n',nfuneval+evals_local,fbest,cputime-cpu_time);
                                end

                                if ftemp2<ftemp
                                    xtemp=xtemp2;
                                    ftemp=ftemp2;
                                else
                                    parar=1;

                                end
                            end
                        end
                        %end of go beyond

                        if search_type==1 %First improvement
                            xcurrloc=xtemp;
                            fcurrloc=ftemp;
                            last_var_used=u_vars(ccc(i));
                            exit_local=1;
                            %                             pause

                            break       %VER SI ESTE BREAK SALE DE LOS DOS FOR
                        else %Best improvement
                            x_improvements=[x_improvements;xtemp];
                            f_improvements=[f_improvements;ftemp];
                            vars_improving=[vars_improving;u_vars(ccc(i))];
                        end
                    end

                    if (nfuneval+evals_local)>=maxevals | cputime-cpu_time>=maxtime
                        improve=improve_solution;
                        if search_type==2 & local_improvement
                            [best_child,iii]=min(f_improvements);
                            xout=x_improvements(iii,:);
                            fout=best_child;
                        else 
                            xout=xcurrloc;
                            fout=fcurrloc;
                        end
                        return
                    end
                end
            end
            if exit_local,  break; end
        end
    end
 
    if ~local_improvement, fin=1; end       %Get out if the local search did not improve the current solution

    if local_improvement & search_type==2 %Once we got out from the 2nd for loop, we check if we explore all the possible u_vars
        [best_child,iii]=min(f_improvements);
        xcurrloc=x_improvements(iii,:);
        fcurrloc=best_child;
        last_var_used=vars_improving(iii);
    end
end

xout=xcurrloc;
fout=fcurrloc;

improve=improve_solution;



