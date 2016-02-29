function print(m,options,pr2lat,pw)
global pwGlobals
    if~exist('options','var')||isempty(options)
        options='org';
    end
    
    if~exist('pr2lat','var')||isempty(pr2lat)
        pr2lat=0;
    end
    
    if~exist('pw','var')||isempty(pw)
       pw=0; 
    end


table=makeTable(m);

%% which option ?

    if regexp(options,'org')
        ix=1:length(m.S(:,1));
        printTable(table,ix);
       
        if(pr2lat==1)
        print2latex(table,ix);
        print2latexOnlyS(table,ix);
        end
    end

    if regexp(options,'r2')
        [s,ix]=sort(table(:,end-2));
        % reorder ix
        ix=ix(end:-1:1);
        
        table=table(ix,:);
        printTable(table,ix);
        if(pr2lat==1)
        print2latex(table,ix);
        print2latexOnlyS(table,ix);
        end
    end
    
    if regexp(options,'cv')
        [s,ix]=sort(table(:,end-1));
        % reorder ix
        ix=ix(end:-1:1);
        
        table=table(ix,:);
        printTable(table,ix);
        if(pr2lat==1)
        print2latex(table,ix);
        print2latexOnlyS(table,ix);
        end
    end
    
    
    if regexp(options,'#')
        [s,ix]=sort(table(:,end));
        % reorder ix
        ix=ix(end:-1:1);
        
        table=table(ix,:);
        printTable(table,ix);
        if(pr2lat==1)
        print2latex(table,ix);
        print2latexOnlyS(table,ix);
        end
    end
    
if(pw==1)
    
IDs = pwGlobals.parsForFitIDs(pwGlobals.indFittedPars);
if(isempty(pwGlobals.indFixedPars)==0)
%cut out from IDs
      IDsStorage=cell(1);
      for i=1:length(IDs)
         if(sum(i==pwGlobals.indFixedPars)~=1) 
            IDsStorage(end+1)=IDs(i);
         end
      end
      IDs=IDsStorage(2:end);
end

for i=1:length(IDs)
    if i<10
        par=(['p' num2str(i) '  ' cell2mat(IDs(i))]); 
    else
        par=(['p' num2str(i) ' ' cell2mat(IDs(i))]); 
    end
    
    disp(sprintf('%s',par))
end

disp(' ')
    
end

end



%% output

function printTable(table,ix)

threshold_r2=0.9;
threshold_cv=0.1;

disp(' ')
%Col-names
    a=sprintf(' p%i ',1);
    for i=2:length(table(:,1))
        if ix(i)<10
            a=[a ' ' sprintf('p%i ',i)];
        else
            a=[a ' ' sprintf('p%i',i)];
        end
        
    end
    a=['ix  |'   a ' |   r2' '    cv' '    #' '  pars'];
    disp(a)
    
% numerical output
for i=1:length(table(:,1))
    
    % index
    if ix(i)<10
        out_ix=sprintf('%i ',ix(i));
    else
        out_ix=sprintf('%i',ix(i));
    end
    
    out_S=sprintf(' %i  ',table(i,1:end-3));
    
    out_r2_cv=sprintf('%1.3f  ',table(i,end-2:end-1));
    
    out_num=sprintf('%i ',table(i,end));
    
    out_pars=sprintf('p%i,',ToPn(table(i,1:end-3),'1'));
    placesOfCommas=regexp(out_pars,',');
    out_pars(placesOfCommas(end))=[];
    
    out_star=[];
    if (table(i,end-2)>threshold_r2) && (table(i,end-1)>threshold_cv) 
        out_star=[out_star '*'];
        if (table(i,end)>1)
            out_star=[out_star '*'];
        end
    end
   
        
        disp([out_ix '  |' out_S ' | ' out_r2_cv '' out_num ' ' out_pars ' ' out_star])

    
end
    disp(' ')
    disp('(*) r2>0.9 & cv > 0.1 (**) r2>0.9 & cv > 0.1 & #>1')
    disp(' ')
end

%% print to latex
function print2latex(table,ix)

threshold_r2=0.9;
threshold_cv=0.1;

fid=fopen('motaPrintOut.txt','a+');
c='c|ccccc';
% for i=2:length(table(1,:))
%    c=[c 'c']; 
% end
fprintf(fid,['\\begin{tabular}{' c '}\n']);

%Col-names
%     a=sprintf('p%i ',1);
%     for i=2:length(table(:,1))
%         if ix(i)<10
%             a=[a ' ' sprintf('&p%i ',i)];
%         else
%             a=[a ' ' sprintf('&p%i',i)];
%         end
%         
%     end
%   fprintf(fid,['ix&'   a ' r2 & cv & $\\sharp$&  pars& \\\\ \\hline' '\n']);
    fprintf(fid,['ix & r2 & cv & $\\sharp$&  pars& scr \\\\ \\hline' '\n']);
    
% numerical output
for i=1:length(table(:,1))
    
    % index
    if ix(i)<10
        out_ix=sprintf('%i ',ix(i));
    else
        out_ix=sprintf('%i',ix(i));
    end
    
%     out_S=sprintf(' %i&  ',table(i,1:end-3));
    
    out_r2_cv=sprintf('%1.3f&  ',table(i,end-2:end-1));
    
    out_num=sprintf('%i ',table(i,end));
    
    out_pars=sprintf('p%i,',ToPn(table(i,1:end-3),'1'));
    placesOfCommas=regexp(out_pars,',');
    out_pars(placesOfCommas(end))=[];
    
    out_star=[];
    if (table(i,end-2)>threshold_r2) && (table(i,end-1)>threshold_cv) 
        out_star=[out_star '*'];
        if (table(i,end)>1)
            out_star=[out_star '*'];
        end
    end
   
        
%       fprintf(fid,[out_ix '&' out_S  out_r2_cv  out_num '&' out_pars '&' out_star '\\\\' '\n']);
        fprintf(fid,[out_ix '&' out_r2_cv  out_num '&' out_pars '&' out_star '\\\\' '\n']);
    
end
fprintf(fid,'\\end{tabular}\\\\');
fclose(fid);
disp('table printed to motaPrintOut.txt')
end

%% 
function print2latexOnlyS(table,ix)
fid=fopen('motaPrintOutOnlyS.txt','a+');
c='c';
for i=2:length(table(1,:))-3
   c=[c 'c']; 
end
fprintf(fid,['\\begin{eqnarray*}\\left( \\begin{array}{' c '}\n']);
%Col-names
    a=sprintf(' p%i ',1);
    for i=2:length(table(:,1))
        if ix(i)<10
            a=[a ' ' sprintf('&p%i ',i)];
        else
            a=[a ' ' sprintf('&p%i',i)];
        end
        
    end
    fprintf(fid,[a ' \\\\ ' '\n']);

for i=1:length(table(:,1))
    out_S=sprintf('&%i  ',table(i,1:end-3));
    out_S=out_S(2:end);
    fprintf(fid,[out_S '\\\\' '\n']);
end

fprintf(fid,'\\end{array}\\right)\\end{eqnarray*}\\\\');
fclose(fid);
disp('matrix S printed to motaPrintOutOnlyS.txt')


end

%% ToPn
% Frequently Used. It converts a vector of ones and zeros to the
% corresponding Number of the Parameter. For Example:
% [0 1 0 1 0 0 0 1] -> [2 4 8]
function Pn=ToPn(IUL,s)
Pn=(regexp(num2str(IUL),s')+2)/3;
end