% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_aux_local.m 1804 2014-07-14 14:32:17Z attila $
function ssm_aux_local(fobj,fjac,x_L,x_U,c_L,c_U,neq,local_solver,nvar,varargin)

params=varargin;

global n_upper n_lower ccuu ccll neqc nvarc

nvarc=nvar;
neqc=neq;

n_upper=[];
n_lower=[];
ccll=[];
ccuu=[];


switch local_solver
    case 'nomad'
        if isempty(c_U)
            n_upper=[];
            n_lower=[];
        else
            %Number of nlc's
            n_upper=[1:length(c_U)];
            n_lower=n_upper;
            clinf=find(c_L==-inf);
            cuinf=find(c_U==inf);

            %Upper bounded constraints
            n_upper(cuinf)=[];

            %Lower bounded contraints constraints
            n_lower(clinf)=[];
        end


        %File with the parameters (for constrained problems)
        if not(isempty(c_U))% | not(isempty(params))
            fid=fopen('fobj_nomad_Param.m','w');
            fprintf(fid,'function Param= f_Param\n');

            fprintf(fid,'global params \n');
            if not(isempty(c_U))
                fprintf(fid,'Param.n_lower =[');
                for i=1:length(n_lower)
                    fprintf(fid,' %g;',n_lower(i));
                end
                fprintf(fid,' ];\n');

                fprintf(fid,'Param.n_upper =[');
                for i=1:length(n_upper)
                    fprintf(fid,' %g;',n_upper(i));
                end
                fprintf(fid,' ];\n');
                fprintf(fid,'Param.c_L =[');
                for i=1:size(c_L,2)
                    fprintf(fid,' %g;',c_L(i));
                end
                fprintf(fid,' ];\n');
                fprintf(fid,'Param.c_U =[');
                for i=1:size(c_U,2)
                    fprintf(fid,' %g;',c_U(i));
                end
                fprintf(fid,' ];\n');
            end

            if not(isempty(params))
                fprintf(fid,'Param.params =params; \n');
                %                 for i=1:length(params)
                %                     fprintf(fid,' %f;',params{i});
                %                 end
%                fprintf(fid,' };\n');
            end
            fclose(fid);
        end

        %File with objective function and constraints
        fid=fopen('fobj_nomad.m','w');
        if isempty(c_U)
            fprintf(fid,'function [fx]= fobj_nomad (x)\n');
        else
            fprintf(fid,'function [fx,cx]= fobj_nomad (x)\n');
            fprintf(fid,'Param=getappdata(0, ''PARAM'' );\n');
        end

        %If there are constraints
        if not(isempty(c_U))
            if not(isempty(params))
                fprintf(fid,'global input_par \n');
                fprintf(fid,'[fx ggg]= %s(x,input_par{:});\n',fobj);
            else
                fprintf(fid,'[fx ggg]= %s(x);\n',fobj);
            end

            fprintf(fid,'cx=[];\n');
            fprintf(fid,'for i=1:length(Param.n_upper) \n');
            fprintf(fid,'cx=[cx ggg(Param.n_upper(i))-Param.c_U(Param.n_upper(i))]; \n');
            fprintf(fid,'end \n');

            fprintf(fid,'for j=1:length(Param.n_lower) \n');
            fprintf(fid,'cx=[cx Param.c_L(Param.n_lower(j))-ggg(Param.n_lower(j))]; \n');
            fprintf(fid,'end \n');
        else
            if not(isempty(params))
                fprintf(fid,'global input_par \n');
                fprintf(fid,'[fx]= %s(x,input_par{:});\n',fobj);
            else
                fprintf(fid,'[fx]= %s(x);\n',fobj);
            end
        end

        fprintf(fid,'return\n');
        fclose(fid);

        %File with bounds
        fid=fopen('fobj_nomad_Omega.m','w');
        fprintf(fid,'function [A,l,u]= fobj_nomad_Omega (n)\n');
        fprintf(fid,'A = [eye(n)];\n');

        fprintf(fid,'l= [');
        for i=1:size(x_L,2)
            fprintf(fid,' %g;',x_L(i));
        end

        fprintf(fid,' ];\n');
        fprintf(fid,'u= [');
        for i=1:size(x_U,2)
            fprintf(fid,' %g;',x_U(i));
        end
        fprintf(fid,' ];\n');

        fprintf(fid,'return\n');
        fclose(fid);



    case {'n2fb','dn2fb'}
        %Creating objective function for n2fb and dn2fb
        if nargout(fobj)<3
            error('The objective function must have at least 3 arguments when using n2fb or dn2fb')
        end
        switch local_solver
            case 'n2fb'
                fid=fopen('objf_n2fb.m','w');
                fprintf(fid,'function [R]= objf_n2fb (N,P,x,NF,R,LTY,TY)\n');
            case 'dn2fb'
                fid=fopen('objf_dn2fb.m','w');
                fprintf(fid,'function [R]= objf_dn2fb (N,P,x,NF,R,LTY,TY)\n');
        end
        fprintf(fid,'global n_fun_eval \n');

        if not(isempty(params))
            fprintf(fid,'global input_par \n');
            fprintf(fid,'[f,g,R]= %s(x,input_par{:});\n',fobj);
        else
            fprintf(fid,'[f,g,R]= %s(x);\n',fobj);
        end

        fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
        fprintf(fid,'return\n');
        fclose(fid);
    case 'nl2sol'
        %Creating objective function for  nl2sol
        if nargout(fobj)<3
            error('The objective function must have at least 3 arguments when using nl2sol: F, constr, R')
        end
        % see if the 
        flag = ess_aux_local_checkfiles('objf_nl2sol.m', fobj);
        if flag == 0
            fid=fopen('objf_nl2sol.m','w');
            fprintf(fid,'function [R]= objf_nl2sol(x)\n');
            fprintf(fid,'%% Objective function for NL2SOL in AMIGO\n%% Automatically generated in ssm_aux_local.m\n');
            fprintf(fid,'global n_fun_eval \n');
            
            if not(isempty(params))
                fprintf(fid,'global input_par \n');
                fprintf(fid,'[f,g,R]= %s(x,input_par{:});\n',fobj);
            else
                fprintf(fid,'[f,g,R]= %s(x);\n',fobj);
            end
            
            fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
            fprintf(fid,'return\n');
            fclose(fid);
        end

        % creating the gradient function if defined
        if ~isempty(fjac)
            flag = ess_aux_local_checkfiles('fjac_nl2sol.m', fjac);
            if flag == 0
                fid=fopen('fjac_nl2sol.m','w');
                fprintf(fid,'function [Jres]= fjac_nl2sol(x)\n');
                fprintf(fid,'%% Jacobian of objective function for NL2SOL in AMIGO\n%% Automatically generated in ssm_aux_local.m\n');
                fprintf(fid,'global n_jac_eval \n');
                
                if not(isempty(params))
                    fprintf(fid,'global input_par \n');
                    fprintf(fid,'[Jobj, Jg, Jres]= %s(x,input_par{:});\n',fjac);
                else
                    fprintf(fid,'[Jobj, Jg, Jres]= %s(x);\n',fjac);
                end
                %fprintf(fid,'[J]= %s(x,inputs);\n',fjac);
                fprintf(fid,'n_jac_eval=n_jac_eval+1; \n\n');
                fprintf(fid,'return\n');
                fclose(fid);
            end
        end

    case 'fsqp'
        %File with the objective function and constraints
        fid=fopen('fobj_fsqp.m','w');
        fprintf(fid,'function [f]= fobj_fsqp (x,j)\n');
        fprintf(fid,'global n_fun_eval \n');
        if not(isempty(params))
            fprintf(fid,'global input_par \n');
            fprintf(fid,'[f]= %s(x,input_par{:});\n',fobj);
        else
            fprintf(fid,'[f]= %s(x);\n',fobj);
        end

        fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
        fprintf(fid,'return\n');
        fclose(fid);

        if not(isempty(c_U))
            ccuu=c_U;
            ccll=c_L;

            if not(isempty(c_U))
                %Number of nlc's
                n_upper=[1:length(c_U)];
                n_lower=n_upper;
                clinf=find(c_L==-inf);
                cuinf=find(c_U==inf);

                %Upper bounded constraints
                n_upper(cuinf)=[];

                %Lower bounded constraints
                n_lower(clinf)=[];
            else
                n_upper=[];
                n_lower=[];
            end




            %File with the objective function and constraints
            fid=fopen('constr_fsqp.m','w');
            fprintf(fid,'function gj= constr_fsqp(x,j)\n');
            fprintf(fid,'global n_upper n_lower ccuu ccll neqc \n');
            fprintf(fid,'global n_fun_eval \n');
            if not(isempty(params))
                fprintf(fid,'global input_par \n');
                fprintf(fid,'[fx ggg]= %s(x,input_par{:});\n',fobj);
            else
                fprintf(fid,'[fx ggg]= %s(x);\n',fobj);
            end

            fprintf(fid,'Amatrix=[];\n');
            fprintf(fid,'c=[];\n');
            fprintf(fid,'for i=1:length(n_upper) \n');
            fprintf(fid,'c=[c ggg(n_upper(i))-ccuu(n_upper(i)-neqc)]; \n');
            fprintf(fid,'end \n');

            fprintf(fid,'for k=1:length(n_lower) \n');
            fprintf(fid,'c=[c ccll(n_lower(k)-neqc)-ggg(n_lower(k))]; \n');
            fprintf(fid,'end \n');



            fprintf(fid,'for i=1:neqc \n');
            fprintf(fid,'c=[c ggg(i)]; \n');
            fprintf(fid,'end \n');

            fprintf(fid,'gj=c(j); \n');
            fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
            %fprintf(fid,'keyboard \n');
            %fprintf(fid,'ggg\n');
            %fprintf(fid,'c \n');
            fprintf(fid,'return\n');
            fclose(fid);
        end

    case 'ipopt'
        %         %File with the objective function and constraints
        %         fid=fopen('ipopt_f.m','w');
        %         fprintf(fid,'function f= ipopt_f(x,Params)\n');
        %         fprintf(fid,'global n_fun_eval \n');
        %         fprintf(fid,'params=Params.params;\n');
        %         fprintf(fid,'f= %s(x,params{:});\n',fobj);
        %         fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
        %         fprintf(fid,'return\n');
        %         fclose(fid);



        %File with the objective function and constraints
        fid=fopen('ipopt_f.m','w');
        fprintf(fid,'function f=ipopt_f(x)\n');
        fprintf(fid,'global n_fun_eval params \n');
        fprintf(fid,'f= %s(x,params{:});\n',fobj);
        fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
        fprintf(fid,'return\n');
        fclose(fid);




        if not(isempty(c_U))
            ccuu=c_U;
            ccll=c_L;
            if isempty(c_U)
                n_upper=[];
                n_lower=[];
            else
                %Number of nlc's
                n_upper=[1:length(c_U)];
                n_lower=n_upper;
                clinf=find(c_L==-inf);
                cuinf=find(c_U==inf);

                %Upper bounded constraints
                n_upper(cuinf)=[];

                %Lower bounded constraints
                n_lower(clinf)=[];
            end
        

        
            fid2=fopen('ipopt_c.m','w');
            fprintf(fid2,'function [c]= ipopt_f(x)\n');
            fprintf(fid,'global n_fun_eval params \n');
            fprintf(fid2,'[f,c]= %s(x,params{:});\n',fobj);
            fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
            fprintf(fid2,'return\n');
            fclose(fid2);
        end


        %             fid2=fopen('ipopt_c.m','w');
        %             fprintf(fid2,'function ceq= ipopt_f(x,dummy)\n');
        %             %We could use input parameters instead of global variables?
        %             fprintf(fid2,'global n_upper n_lower ccuu ccll nvarc neqc \n');
        %             fprintf(fid,'global n_fun_eval \n');
        %             fprintf(fid,'params=dummy.params;\n');
        %             fprintf(fid2,'[f,ggg]= %s(x,params{:});\n',fobj);
        %             fprintf(fid2,'ceq=[];\n');
        %             fprintf(fid2,'ceq=[ceq ggg(1:neqc)];\n');
        %             fprintf(fid2,'clow=[];\n');
        %             fprintf(fid2,'cupp=[];\n');
        %             fprintf(fid2,'for i=1:length(n_upper) \n');
        %             fprintf(fid2,'cupp=[cupp; ggg(n_upper(i))-ccuu(n_upper(i)-neqc)+x(nvarc+i)]; \n');
        %             fprintf(fid2,'end \n');
        %             fprintf(fid2,'for j=1:length(n_lower) \n');
        %             fprintf(fid2,'clow=[clow; ggg(n_lower(j))-ccll(n_lower(j)-neqc)-x(nvarc+length(n_upper)+j)]; \n');
        %             fprintf(fid2,'end \n');
        %             % fprintf(fid2,'keyboard \n');
        %             fprintf(fid2,'ceq=[ceq''; cupp; clow]; \n');
        %
        %             fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
        %             fprintf(fid2,'return\n');
        %             fclose(fid2);



case 'misqp'
    if not(isempty(c_U))
        ccuu=c_U;
        ccll=c_L;
        if isempty(c_U)
            n_upper=[];
            n_lower=[];
        else
            %Number of nlc's
            n_upper=[1:length(c_U)];
            n_lower=n_upper;
            clinf=find(c_L==-inf);
            cuinf=find(c_U==inf);

            %Upper bounded constraints
            n_upper(cuinf)=[];

            %Lower bounded constraints
            n_lower(clinf)=[];
        end
    end
        case 'lbfgsb'
        %Creating objective function for  lbfgsb
        if nargout(fobj)< 3
            error('The objective function must have 3 output arguments (f_obj 0 0)')
        end
        
        fid=fopen('objf_lbfgsb.m','w');
        fprintf(fid,'function [f]= objf_lbfgsb(x)\n');

        fprintf(fid,'global n_fun_eval \n');

        if not(isempty(params))
            fprintf(fid,'global input_par \n');
            fprintf(fid,'[f,g,R]= %s(x,input_par{:});\n',fobj);
        else
            fprintf(fid,'[f,g,R]= %s(x);\n',fobj);
        end

        fprintf(fid,'n_fun_eval=n_fun_eval+1; \n');
        fprintf(fid,'return\n');
        fclose(fid);
        
        
        % creating the gradient function if defined
        if ~isempty(fjac)
            fid=fopen('fjac_lbfgsb.m','w');
            fprintf(fid,'function [Jobj]= fjac_lbfgsb(x)\n');
            fprintf(fid,'%% Automatically generated in ssm_aux_local.m\n');
            fprintf(fid,'global n_jac_eval \n');
            fprintf(fid,'global input_par \n\n');
            fprintf(fid,'[Jobj, Jg, Jres]= %s(x,input_par{:});\n',fjac);           
            fprintf(fid,'n_jac_eval=n_jac_eval+1; \n\n');
            fprintf(fid,'return\n');
            fclose(fid);
        end

    otherwise
        if not(isempty(c_U))
            ccuu=c_U;
            ccll=c_L;
            if isempty(c_U)
                n_upper=[];
                n_lower=[];
            else
                %Number of nlc's
                n_upper=[1:length(c_U)];
                n_lower=n_upper;
                clinf=find(c_L==-inf);
                cuinf=find(c_U==inf);

                %Upper bounded constraints
                n_upper(cuinf)=[];

                %Lower bounded constraints
                n_lower(clinf)=[];
            end

        end
end
rehash

n_upper(1:neqc)=[];
n_lower(1:neqc)=[];
%ccll(1:neqc)=[];
%ccuu(1:neqc)=[];
