function [transf_new,nv_new] = sm(transf,nuevas_variables,allVar)
    assume(transf,'real')
    ttransf=transf;
    [m,n]=size(ttransf);
    ind=[];
    for i=1:m
       for j=i+1:m
          aux=dot(transf(i,:),transf(j,:));
          if aux==0
              ind=[ind;i,j];
          end
       end
    end
    [m1,n1]=size(ind);

    fprintf('Generadores:\n')
    disp(transf)
    if m1~=0
        transf_new=transf;
        nv_new=nuevas_variables;
        for i=1:m1
           disp(transf(ind(i,1),:)+transf(ind(i,2),:))
           transf_new=[transf_new;transf(ind(i,1),:)+transf(ind(i,2),:)];
%            k=find(transf_new(end,:));
           AV=allVar';
           aux=subs(AV,allVar',nuevas_variables(ind(i,1),:));
           aux=subs(aux,allVar',nuevas_variables(ind(i,2),:));
           nv_new=[nv_new;aux];
        end
        
    else
        transf_new=transf;
        nv_new=nuevas_variables;
    end
    
 
end

