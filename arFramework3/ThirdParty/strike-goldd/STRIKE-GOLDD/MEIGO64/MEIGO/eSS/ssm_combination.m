% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_combination.m 770 2013-08-06 09:41:45Z attila $
function [new_comb,new_comb_values, new_comb_values_penalty,new_comb_nlc, new_comb_penalty,...
    nfuneval,n_new]=ssm_combination(x1,x2,...
    x1_val,x2_val,nrand,x_L,x_U,nconst,fobj,c_L,c_U,tolc,weight,nfuneval,int_var,bin_var,dim_refset,index,prob_bound,varargin);

nvar=length(x1);
%Maximum number of solutions evaluated = 100 by now
new_comb=zeros(100,nvar);
new_comb_values=zeros(100,1);
new_comb_values_penalty=zeros(100,1);
new_comb_penalty=zeros(100,1);
new_comb_nlc=zeros(100,max(1,nconst));

n_new=1;


c=[];
p=[];
pval=[];


%Auxliary vectors
v=zeros(5,nvar);
d=(x2-x1)/2;
v(1,:)=x1-d;

aaa=find(v(1,:)<x_L);
bbb=find(v(1,:)>x_U);

if ~isempty(aaa)
    if rand>prob_bound
        v(1,aaa)=x_L(aaa);
    end
end

if ~isempty(bbb)
    if rand>prob_bound
        v(1,bbb)=x_U(bbb);
    end
end

v(2,:)=x1;

v(3,:)=x1+d;

aaa=find(v(3,:)<x_L);
bbb=find(v(3,:)>x_U);

if ~isempty(aaa)
    if rand>prob_bound
        v(3,aaa)=x_L(aaa);
    end
end

if ~isempty(bbb)
    if rand>prob_bound
        v(3,bbb)=x_U(bbb);
    end
end


v(4,:)=x2;
v(5,:)=x2+d;

aaa=find(v(5,:)<x_L);
bbb=find(v(5,:)>x_U);

if ~isempty(aaa)
    if rand>prob_bound
        v(5,aaa)=x_L(aaa);
    end
end

if ~isempty(bbb)
    if rand>prob_bound
        v(5,bbb)=x_U(bbb);
    end
end



%If both elements belong to Refset1
if (index(1)<=dim_refset/2) & (index(2)<=dim_refset/2)
    c(1,:)=v(1,:)+(v(2,:)-v(1,:)).*rand(1,nrand);   %C1
    p(1,:)=x1;
    pval(1)=x1_val;
    c(2,:)=v(2,:)+(v(3,:)-v(2,:)).*rand(1,nrand);   %C2
    p(2,:)=x1;
    pval(2)=x1_val;
    c(3,:)=v(4,:)+(v(5,:)-v(4,:)).*rand(1,nrand);   %C3
    p(3,:)=x2;
    pval(3)=x2_val;

    %2nd C2
    c(4,:)=v(3,:)+(v(4,:)-v(3,:)).*rand(1,nrand);
    p(4,:)=x2;
    pval(4)=x2_val;


    n_combin=4;

    %%If both elements belong to Refset2
elseif (index(1)>dim_refset/2) & (index(2)>dim_refset/2)
    c(1,:)=v(2,:)+(v(3,:)-v(2,:)).*rand(1,nrand);
    p(1,:)=x1;
    pval(1)=x1_val;
    a=rand;
    if a<0.5
        c(2,:)=v(2,:)+(v(1,:)-v(2,:)).*rand(1,nrand);          %C1
        p(2,:)=x1;
        pval(2)=x1_val;
    else
        c(2,:)=v(4,:)+(v(5,:)-v(4,:)).*rand(1,nrand);
        p(2,:)=x2;
        pval(2)=x2_val;
    end
    n_combin=2;

    %One from Refset1, another one from Refset2
else
    c(1,:)=v(2,:)+(v(1,:)-v(2,:)).*rand(1,nrand);   %C1   %C1
    p(1,:)=x1;
    pval(1)=x1_val;
    c(2,:)=v(2,:)+(v(3,:)-v(2,:)).*rand(1,nrand);   %C2
    p(2,:)=x1;
    pval(2)=x1_val;
    c(3,:)=v(4,:)+(v(5,:)-v(4,:)).*rand(1,nrand);   %C3
    p(3,:)=x2;
    pval(3)=x2_val;
    n_combin=3;
end



for i=1:n_combin;

    %Evaluate
    [val,val_penalty,pena,nlc,include,x] = ssm_evalfc(c(i,:),x_L,x_U,fobj,nconst,...
        c_L,c_U,tolc,weight,int_var,bin_var,varargin{:});

    nfuneval=nfuneval+1;


    if include
        new_comb(n_new,:)=x;
        new_comb_values(n_new)=val;
        new_comb_values_penalty(n_new)=val_penalty;
        new_comb_nlc(n_new,:)=nlc;
        new_comb_penalty(n_new)=pena;
        n_new=n_new+1;


        if val_penalty<pval(i) & index(1)<=dim_refset/2
            if not(index(2)>dim_refset/2 & i==3)

                [vector,vector_value,vector_penalty,vector_value_penalty,vector_nlc,new_child,new_child_value,...
                    new_child_penalty,new_child_value_penalty,new_child_nlc,nfuneval]=ssm_beyond(p(i,:),...
                    x,val_penalty,fobj,nrand,tolc,weight,...
                    x_L,x_U,c_L,c_U,nconst,int_var,bin_var,nfuneval,prob_bound,varargin{:});


                new_child_index=length(new_child_value);


                if ~isempty(new_child)
                    new_comb(n_new:n_new+new_child_index-1,:)=new_child;
                    new_comb_values(n_new:n_new+new_child_index-1)=new_child_value;
                    new_comb_values_penalty(n_new:n_new+new_child_index-1)=new_child_value_penalty;
                    new_comb_nlc(n_new:n_new+new_child_index-1,:)=new_child_nlc;
                    new_comb_penalty(n_new:n_new+new_child_index-1)=new_child_penalty;
                    n_new=n_new+new_child_index;
                end
            end
        end
    end
end