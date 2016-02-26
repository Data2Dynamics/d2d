function table=makeTable(m)    
table=NaN(length(m.S(:,1)),length(m.S(:,1))+3);
    table(:,end-1)=std(m.K)./abs(mean(m.K));

    S_str=mat2str(m.S);
    
    for i=1:length(m.S(:,1))
        
        help=find(m.r2(2,:,i));
        if isempty(help)==1
            help=1;
        end
        table(i,end-2)=m.r2(2,help(end),i);
        
        % find number of identical funtional relations
        a=num2str(m.S(i,1));
        for j=2:length(m.S(:,1))
           
            a=[a ' ' num2str(m.S(i,j))]; 
        end
        
        S_regexp=regexp(S_str,a);
        table(i,end)=length(S_regexp);
        
        
    end
    table(:,1:end-3)=m.S;

end
    