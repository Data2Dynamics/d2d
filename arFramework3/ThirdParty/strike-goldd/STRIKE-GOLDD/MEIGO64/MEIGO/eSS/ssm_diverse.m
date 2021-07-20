% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_diverse.m 770 2013-08-06 09:41:45Z attila $
function fx=diverse(ndiverse,x_L,x_U)

n=length(x_L);

%Generates ndiverse diverse vectors between the bounds x_L and x_U

f=ones(n,4);                     %Matrix of frequencies

solutions=zeros(ndiverse,n);          %Matrix of solutions
p=zeros(4,1);                    %Vector of probabilities



%We generate 4 solutions
for i=1:4
solutions(i,:)=(rand(1,n)+i-1)/4;
end

for i=5:ndiverse                           %ndiverse-1 since we will add the user given solution later
    for j=1:n
        for k=1:4
           p(k)=(1/f(j,k))/sum(1./f(j,:));
        end
        a=rand;
        for m=1:4
            if a<=sum(p(1:m))
                solutions(i,j)=(rand+m-1)/4;
                f(j,m)=f(j,m)+1;
                break
            end
        end
    end
end

%Put the variables inside the bounds
a=repmat(x_U-x_L,[ndiverse,1]);
b=repmat(x_L,[ndiverse,1]);


solutions=solutions.*a+b;

fx=solutions;