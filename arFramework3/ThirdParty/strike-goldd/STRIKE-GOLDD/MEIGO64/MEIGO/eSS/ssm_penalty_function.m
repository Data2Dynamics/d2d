% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_penalty_function.m 770 2013-08-06 09:41:45Z attila $
function fx=penalty_function(x,nlc,c_L,c_U,tolc);

P=[];

if ~(isempty(nlc))


    if size(nlc)~=size(c_U), nlc=nlc'; end
    
    a=find(nlc<c_L);
    b=find(nlc>c_U);
    
    P=zeros(1,length(a)+length(b));
    counter=1;
    for i=1:size(a,2)
       P(counter)=(c_L(a(i))-nlc(a(i)));
        %P=[P (c_L(a(i))-nlc(a(i)))];%/(1+abs(c_L(a(i))))];
        counter=counter+1;
   end
   
    for i=1:size(b,2)
        P(counter)=(nlc(b(i))-c_U(b(i)));
   %     P=[P (nlc(b(i))-c_U(b(i)))];%/(1+abs(c_U(b(i))))];
       counter=counter+1;
    end

end

%maxP=max(P);

if ~(isempty(P)) & max(P)>tolc
    %The penalty is the maximum constraint violation
    fx=max(P);
else
    fx=0;
end
return
