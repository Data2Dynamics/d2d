function tex = MOTAResultToTex(m)
%
% Returns a cell array.

tex = {};

threshold_r2=0.9;
threshold_cv=0.1;

table=makeTable(m);
c='c|ccccc';
ix=1:length(m.S(:,1));

%% numerical output

tex{end+1} = ['\begin{tabular}{' c '}'];
tex{end+1} = 'ix & r2 & cv & $\sharp$&  pars& scr \\ \hline';
    
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
        tex{end+1} = [out_ix '&' out_r2_cv  out_num '&' out_pars '&' out_star '\\\\'];
    
end
tex{end+1} = '\end{tabular}\\';
tex{end+1} = '';

%% Matrix S

c='c';
for i=2:length(table(1,:))-3
   c=[c 'c']; 
end

tex{end+1} = ['\begin{eqnarray*}\left( \begin{array}{' c '}'];

%Col-names
    a=sprintf(' p%i ',1);
    for i=2:length(table(:,1))
        if ix(i)<10
            a=[a ' ' sprintf('&p%i ',i)];
        else
            a=[a ' ' sprintf('&p%i',i)];
        end
        
    end
    tex{end+1} = [a ' \\ '];

for i=1:length(table(:,1))
    out_S=sprintf('&%i  ',table(i,1:end-3));
    out_S=out_S(2:end);
    tex{end+1} = [out_S '\\'];
end

tex{end+1} = '\end{array}\right)\end{eqnarray*}\\';

end

%% ToPn
% Frequently Used. It converts a vector of ones and zeros to the
% corresponding Number of the Parameter. For Example:
% [0 1 0 1 0 0 0 1] -> [2 4 8]
function Pn=ToPn(IUL,s)
Pn=(regexp(num2str(IUL),s')+2)/3;
end