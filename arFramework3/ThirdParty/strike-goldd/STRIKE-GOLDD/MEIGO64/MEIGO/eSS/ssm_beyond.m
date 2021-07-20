% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_beyond.m 770 2013-08-06 09:41:45Z attila $
function  [vector,vector_value,vector_penalty,vector_value_penalty,vector_nlc,...
    new_child,new_child_value,new_child_penalty,new_child_value_penalty,new_child_nlc,nfuneval]=ssm_trascender(z1,z2,z2_val,fobj,nrand,tolc,weight,...
    x_L,x_U,c_L,c_U,nconst,int_var,bin_var,nfuneval,prob_bound,varargin);

continuar=1;
denom=1;
n_improve=1;


vector=[];
vector_value=[];
vector_value_penalty=[];
vector_penalty=[];
vector_nlc=[];

new_child=[];
new_child_value=[];
new_child_value_penalty=[];
new_child_nlc=[];
new_child_penalty=[];


while continuar

    %Continue
    d=(z2-z1)/denom;
    zv(1,:)=z2;
    zv(2,:)=z2+d;

    aaa=find(zv(2,:)<x_L);
    bbb=find(zv(2,:)>x_U);

    if ~isempty(aaa)
        if rand>prob_bound
            zv(2,aaa)=x_L(aaa);
        end
    end

    if ~isempty(bbb)
        if rand>prob_bound
            zv(2,bbb)=x_U(bbb);
        end
    end

    xnew=zv(1,:)+(zv(2,:)-zv(1,:)).*rand(1,nrand);   %non convex

    %Evaluate
    [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(xnew,x_L,x_U,fobj,nconst,c_L,c_U,...
        tolc,weight,int_var,bin_var,varargin{:});
    nfuneval=nfuneval+1;
    if include
        new_child=[new_child;x];
        new_child_value=[new_child_value;val];
        new_child_value_penalty=[new_child_value_penalty;val_penalty];
        new_child_nlc=[new_child_nlc;nlc];
        new_child_penalty=[new_child_penalty;pena];

        if val_penalty<z2_val
            z1=z2;
            z2=xnew;
            z2_val=val_penalty;

            vector=x;
            vector_value=val;
            vector_value_penalty=val_penalty;
            vector_penalty=pena;
            vector_nlc=nlc;
            n_improve=n_improve+1;

            if n_improve==2;
                denom=denom/2;
                n_improve=0;
            end
        else
            %If it does not improve, break the loop
            break
        end

    else
        continuar=0;
    end
end
