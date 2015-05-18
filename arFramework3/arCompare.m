% [same,d1,d2] = arCompare(ar1,ar2,pattern,opt)
% 
%   Comparison of two ar structs/objects
% 
%   ar1       an ar struct, e.g. the global variable "ar"
%   ar2       another ar struct, e.g. "ar" from an workspace
%   pattern   string or cell of strings  
%             patterns are evaluated with regexp.m
%             if pattern is provided, then only a subset of cells are
%             evaluated. Providing patterns is reasonable because the
%             structs often only differ by ar.fevals, ar.info, ...
%   silent    [false]
% 
%   d1,d2   the differences, the structures are similar to 'ar' but not equal
%           since the fields can have different type and should indicate
%           the kind of difference.   
% 
% Examples:
% ar2 = ar;
% ar2.qLog10(:)=~ar.qLog10;
% [same,d1,d2]=arCompare(ar,ar2,'qLog10')
% 
% Check only specific fields:
% ar2 = ar;
% ar2.model(1).data(1).qFit = round(rand(size(ar.model(1).data(1).qFit));
% [same,d1,d2]=arCompare(ar,ar2,{'qLog10','qFit','logfitting'})


function [same,d1,d2] = arCompare(ar1,ar2,pattern,silent)
if(~exist('pattern','var') || isempty(pattern))
    pattern = [];
end
if(~exist('silent','var') || isempty(silent))
    silent = false;
end

% [same,d1,d2] = structcmp(ar1,ar2);
[same,d1,d2] = structcmp2(ar1,ar2,pattern);

if ~silent
    printDifference(d1,d2);
end


function printDifference(d1,d2)
fs = allFields(d1,1);
fs = regexprep(fs,'(\w)\.(\w)','$1(1).$2');
fs = regexprep(fs,'(\w)$','$1(1)');
fs = strrep(fs,'(','{');
fs = strrep(fs,')','}');

for i=1:length(fs)
    fsname = ['ar.',fs{i}];
    fsname = strrep(fsname,'{','(');
    fsname = strrep(fsname,'}',')');
    
    eval(['out1 = d1.',fs{i},';']);
    eval(['out2 = d2.',fs{i},';']);
    if(ischar(out1))
        if(~strcmp(out1,out2))
            fprintf('%40s different: %s vs. %s \n',fsname,out1,out2);
        else
            fprintf('%40s different: %s \n',fsname,out1);
        end
    elseif(isnumeric(out1))
        fprintf('%40s different: %f vs. %f \n',fsname,out1,out2);
    end
end


function [same,d1,d2] = structcmp2(s1,s2,pattern)

same = true;

fn1 = allFields(s1,1);
fn2 = allFields(s2,1);

if(~isempty(pattern))
    if(~iscell(pattern))
        pattern = {pattern};
    end
    use1 = zeros(size(fn1));
    use2 = zeros(size(fn2));
    for i=1:length(pattern)
        use1 = use1 | ~cellfun(@isempty,regexp(fn1,pattern{i}));
        use2 = use2 | ~cellfun(@isempty,regexp(fn2,pattern{i}));
    end
    fn1 = fn1(use1);
    fn2 = fn2(use2);
end

Fn1 = regexprep(fn1,'(\w)\.(\w)','$1(1).$2');
Fn2 = regexprep(fn2,'(\w)\.(\w)','$1(1).$2');
Fn1 = regexprep(Fn1,'(\w)$','$1(1)');
Fn2 = regexprep(Fn2,'(\w)$','$1(1)');

Fn1 = strrep(Fn1,'(','{');
Fn1 = strrep(Fn1,')','}');
Fn2 = strrep(Fn2,'(','{');
Fn2 = strrep(Fn2,')','}');


d1 = struct;
d2 = struct;

notInBoth = setxor(Fn1,Fn2);
if(~isempty(notInBoth))
    same = false;
    notIn1 = setdiff(Fn2,Fn1);
    for i=1:length(notIn1)
        eval(['d2.',notIn1{i},' = ''field only in the 2nd struct'';']);
        eval(['d1.',notIn1{i},' = '''';']);
%         d1.(notIn1{i}) = '';
    end
    notIn2 = setdiff(Fn1,Fn2);
    for i=1:length(notIn2)
%         d1.(notIn2{i}) = 'field only in one 1st struct';
%         d2.(notIn2{i}) = '';
        eval(['d1.',notIn2{i},' = ''field only in the 1st struct'';']);
        eval(['d2.',notIn2{i},' = '''';']);
    end
end


fn = intersect(fn1,fn2);
% fn
% if(~isempty(pattern))
%     fn = intersect(fn,pattern);
% end
Fn = regexprep(fn,'(\w)\.(\w)','$1(1).$2');
Fn = regexprep(Fn,'(\w)$','$1(1)');
Fn = strrep(Fn,'(','{');
Fn = strrep(Fn,')','}');

for i=1:length(fn)
    eval(['val1 = s1.',fn{i},';']);
    eval(['val2 = s2.',fn{i},';']);
    
    if(strcmp(class(val1),class(val2))~=1)
        same = false;
        eval(['d1.',Fn{i},' = sprintf(''class = %s'',class(val1));']);
        eval(['d2.',Fn{i},' = sprintf(''class = %s'',class(val2));']);
        
    elseif(isempty(val1) ~= isempty(val2))
        same = false;
        tmp = sprintf('%i x ',size(val1));
        eval(['d1.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
        tmp = sprintf('%i x ',size(val2));
        eval(['d2.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
        
    elseif(length(size(val1)) ~= length( size(val2)))
        same = false;
        tmp = sprintf('%i x ',size(val1));
        eval(['d1.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
        tmp = sprintf('%i x ',size(val2));
        eval(['d2.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
        
    elseif(sum(abs(size(val1) - size(val2) ))~=0)
        same = false;
        tmp = sprintf('%i x ',size(val1));
        eval(['d1.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
        tmp = sprintf('%i x ',size(val2));
        eval(['d2.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
        
    elseif(isnumeric(val1) || islogical(val1))
        if(sum(abs(val1(:)-val2(:)))>0)
            if(length(val1)<=5)
                eval(['d1.',Fn{i},' = sprintf(''%f '',val1);']);
                eval(['d2.',Fn{i},' = sprintf(''%f '',val2);']);
            else
                eval(['d1.',Fn{i},' = sprintf(''arrays different'');']);
                eval(['d2.',Fn{i},' = sprintf(''arrays different'');']);
            end
            same = false;
        end
        if(sum(abs(isnan(val1(:) ) - isnan( val2(:) ) ))>0)
            same = false;
            eval(['d1.',Fn{i},' = sprintf(''%i NaNs'',sum(isnan(val1(:))));']);
            eval(['d2.',Fn{i},' = sprintf(''%i NaNs'',sum(isnan(val2(:))));']);
        end
    elseif(ischar(val1))
        gleich = strcmp(val1,val2);
        if(~gleich)
            same = false;
            val1 = val1(1:min(length(val1),100));
            val2 = val2(1:min(length(val2),100));
            eval(['d1.',Fn{i},' = val1;']);
            eval(['d2.',Fn{i},' = val2;']);
        end
    elseif(iscell(val1))
        if(isempty(val1) && isempty(val2))
        else
            try
                [csame,reason1,reason2] = cellcmp(val1,val2);
            catch
                csame = 1;
                val1
                val2
            end
            if(~csame)
                same = false;
                eval(['d1.',Fn{i},' = reason1;']);
                eval(['d2.',Fn{i},' = reason2;']);
            end
        end
    elseif(isstruct(val1))
        error('case should not occur.')
%         [same2,dd1,dd2] = structcmp(s1.(fn{i}),val2);
%         if(~same2)
%             same = false;
%             eval(['d1.',Fn{i},' = dd1;
%             d2.(fn{i}) = dd2;
%         end
        
    end
end




function [same,d1,d2] = structcmp(s1,s2)
lmax = max(length(s1),length(s2));

if lmax>1
    same = NaN(size(s1));
    d1 = cell(size(s1));
    d2 = cell(size(s2));
    for i=1:min(10,lmax)
        [same(i),d1{i},d2{i}] = structcmp(s1(i),s2(i));
    end
    for i=(lmax+1):length(s1) % only in s1
        same(i) = false;
        d1{i} = s1(i);
    end
    for i=(lmax+1):length(s2)
        same(i) = false;
        d2{i} = s2(i);
    end
    same = all(same);
    
else
    fn1 = fieldnames(s1);
    fn2 = fieldnames(s2);
    
    same = true;
    
    d1 = struct;
    d2 = struct;
    
    notInBoth = setxor(fn1,fn2);
    if(~isempty(notInBoth))
        same = false;
        notIn1 = setdiff(fn2,fn1);
        for i=1:length(notIn1)
            d2.(notIn1{i}) = 'field only in one 2nd struct';
            d1.(notIn1{i}) = '';
        end
        notIn2 = setdiff(fn1,fn2);
        for i=1:length(notIn2)
            d1.(notIn2{i}) = 'field only in one 1st struct';
            d2.(notIn2{i}) = '';
        end
    end
    
    fn = intersect(fn1,fn2);
    
    for i=1:length(fn)
        if(strcmp(class(s1.(fn{i})),class(s2.(fn{i})))~=1)
            same = false;
            d1.(fn{i}) = sprintf('class = %s',class(s1.(fn{i})));
            d2.(fn{i}) = sprintf('class = %s',class(s2.(fn{i})));
            
        elseif(isempty(s1.(fn{i})) ~= isempty(s2.(fn{i})))
            same = false;
            tmp = sprintf('%i x ',size(s1.(fn{i})));
            d1.(fn{i}) = ['size = ',tmp(1:end-2)];
            tmp = sprintf('%i x ',size(s2.(fn{i})));
            d2.(fn{i}) = ['size = ',tmp(1:end-2)];
            
        elseif(length(size(s1.(fn{i}))) ~= length( size(s2.(fn{i}))))
            same = false;
            tmp = sprintf('%i x ',size(s1.(fn{i})));
            d1.(fn{i}) = ['size = ',tmp(1:end-2)];
            tmp = sprintf('%i x ',size(s2.(fn{i})));
            d2.(fn{i}) = ['size = ',tmp(1:end-2)];
            
        elseif(sum(abs(size(s1.(fn{i})) - size(s2.(fn{i})) ))~=0)
            same = false;
            tmp = sprintf('%i x ',size(s1.(fn{i})));
            d1.(fn{i}) = ['size = ',tmp(1:end-2)];
            tmp = sprintf('%i x ',size(s2.(fn{i})));
            d2.(fn{i}) = ['size = ',tmp(1:end-2)];
            
        elseif(isnumeric(s1.(fn{i})))
            if(sum(abs(s1.(fn{i})(:)-s2.(fn{i})(:)))>0)
                if(length(s1.(fn{i}))<=5)
                    d1.(fn{i}) = sprintf('%f ',s1.(fn{i}));
                    d2.(fn{i}) = sprintf('%f ',s2.(fn{i}));
                else
                    d1.(fn{i}) = sprintf('arrays different');
                    d2.(fn{i}) = sprintf('arrays different');
                end
                same = false;
            end
            if(sum(abs(isnan( s1.(fn{i})(:) ) - isnan( s2.(fn{i})(:) ) ))>0)
                same = false;
                d1.(fn{i}) = sprintf('%i NaNs',sum(isnan(s1.(fn{i})(:))));
                d2.(fn{i}) = sprintf('%i NaNs',sum(isnan(s2.(fn{i})(:))));
            end
        elseif(ischar(s1.(fn{i})))
            gleich = strcmp(s1.(fn{i}),s2.(fn{i}));
            if(~gleich)
                same = false;
                d1.(fn{i}) = s1.(fn{i});
                d2.(fn{i}) = s2.(fn{i});
            end
        elseif(iscell(s1.(fn{i})))
            if(isempty(s1.(fn{i})) && isempty(s2.(fn{i})))
            else
                [csame,reason1,reason2] = cellcmp(s1.(fn{i}),s2.(fn{i}));
                if(~csame)
                    same = false;
                    d1.(fn{i}) = reason1;
                    d2.(fn{i}) = reason2;
                end
            end
        elseif(isstruct(s1.(fn{i})))
            [same2,dd1,dd2] = structcmp(s1.(fn{i}),s2.(fn{i}));
            if(~same2)
                same = false;
                d1.(fn{i}) = dd1;
                d2.(fn{i}) = dd2;
            end
            
        end
    end
    
end


% same = cellcmp(c1,c2)
%
%   Vergleich zwei cell Objekte.
function [same,reason1,reason2] = cellcmp(c1,c2)
same = true;
reason1 = '';
reason2 = '';

if ~isempty(c1)
    isc = find(~cellfun(@isnumeric,c1));
elseif isempty(c1) & isempty(c2)
    return
else
    same = false;
    reason1 = sprintf('isempty=',isempty(c1));
    reason2 = sprintf('isempty=',isempty(c2));
    return
end

if(sum(abs(size(c1)-size(c2)))~=0)
    same = false;
    reason1 = sprintf('size=%i,',size(c1));
    reason2 = sprintf('size=%i,',size(c2));
    return;
elseif sum(abs(cellfun(@iscell,c1(isc))-cellfun(@iscell,c2(isc))))>0
    reason1 = sprintf('cellfun(@iscell,...) different');
    reason2 = sprintf('cellfun(@iscell,...) different');
    same = false;
elseif size(c1,1)==1 && sum(strcmpi(cellfun(@class,c1,'UniformOutput',false), cellfun(@class,c2,'UniformOutput',false))~=1)>0
    reason1 = sprintf('cellfun(@class,...) different');
    reason2 = sprintf('cellfun(@class,...) different');
    same = false;
elseif size(c1,2)==1 && sum(strcmpi(cellfun(@class,c1,'UniformOutput',false)', cellfun(@class,c2,'UniformOutput',false)')~=1)>0
    reason1 = sprintf('cellfun(@class,...) different');
    reason2 = sprintf('cellfun(@class,...) different');
    same = false;
elseif sum(cellfun(@isempty,c1(isc)) - cellfun(@isempty,c2(isc)))>0
    reason1 = sprintf('cellfun(@isempty,...) different');
    reason2 = sprintf('cellfun(@isempty,...) different');
    same = false;
elseif(sum(abs(celllength(c1(isc))-celllength(c2(isc))))>0)
    reason1 = sprintf('length of cell entries');
    reason2 = sprintf('length of cell entries');
    same = false;
    return;
elseif(isempty(c1) && isempty(c2))
elseif(isnumeric(c1{1}))
    c1 = cell2array(c1);
    c2 = cell2array(c2);
    if(sum(abs(c1(:)-c2(:)))~=0)
        reason1 = sprintf('array entries');
        reason2 = sprintf('array entries');
        same = false;
        return;
    end
elseif(ischar(c1{1}))
    [same] = strcmpCell(c1,c2,1);
    reason1 = 'cell of strings different';
    reason2 = 'cell of strings different';
elseif(iscell(c1{1}))
    for i=1:length(c1)
        same = cellcmp(c1{i},c2{i});
        if(~same)
            reason1 = 'nested cell different';
            reason2 = 'nested cell different';
            return;
        end
    end
elseif isstruct(c1)
    same = structcmp(c1,c2);
else
    save error
    class(c1)
    class(c2)
    error('cellcmp(c1,c2) not yet implemented for this class.')
end


% match = strcmpCell(str1,str2)
%
% Vergleicht str1{i} mit str2{j}
%
%   str1        Zelle von Strings z.B. str1 = {'irgendeinstring'}
%   str2        Zelle von Strings z.B. str1 = {'irgendeinstring'}
%   onlyDiagonal    if 1, dann wird nur str1{i} mit str2{i} verglichen
%
%   match       ist eine Matrix length(str1) x length(str2) mit Nullen und Einsen.

function match = strcmpCell(str1,str2,onlyDiagonal)
if(~exist('onlyDiagonal','var') || onlyDiagonal~=1)
    match = NaN(100);
    for i=1:length(str1)
        for j=1:length(str2)
            match(i,j) = strcmp(str1{i},str2{j});
        end
    end
    match = match(1:i,1:j);
else
    match = NaN(1,1000);
    for i=1:length(str1)
        match(i) = strcmp(str1{i},str2{i});
    end
    match = match(1:i);
end



% a = cell2array(c)
%
%     KOnvertiert eine Zelle in einen Array gleicher Dimension.
function a = cell2array(c)
a = NaN*ones(size(c));
c = {c{:}};
cl = celllength(c);

if(sum(cl==1)==length(cl))
    for i=1:length(c)
        try
            a(i) = c{i};
        catch
            c{i}
            length(c{i})
            class(c{i})
        end
    end
elseif(sum((cl==1) | cl==0)==length(cl))
    for i=1:length(c)
        if(~isempty(c{i}))
            a(i) = c{i};
        else
            a(i) = NaN;
        end
    end
else
    a = NaN*ones(1,sum(cl));
    ind = 0;
    for i=1:length(c)
        a(ind+(1:cl(i))) = c{i};
        ind = ind+cl(i);
    end
end


function fns = allFields(s,onlyFinal)
if(~exist('onlyFinal','var') || isempty(onlyFinal))
    onlyFinal = 0;
end

fns = cell(50000,1);
ii=0;

for is = 1:length(s)
    if(isstruct(s(is)))
        
        fn = fieldnames(s(is));
        for i=1:length(fn)
            if(isstruct(s(is).(fn{i})))
                for j=1:length(s(is).(fn{i}))
                    fns_tmp = allFields(s(is).(fn{i})(j));
                    for k=1:length(fns_tmp)
                        ii=ii+1;
                        fns{ii} = [fn{i},'(',int2str(j),').',fns_tmp{k}];
%                         fns{ii} = sprintf('%s(%i).%s',fn{i},j,fns_tmp{k});
                    end
                end
            elseif(iscell(s(is).(fn{i})) && ~isempty(s(is).(fn{i})) && isstruct(s(is).(fn{i}){1}))
                for j=1:length(s(is).(fn{i}))
                    fns_tmp = allFields(s(is).(fn{i}){j});
                    for k=1:length(fns_tmp)
                        ii=ii+1;
                        fns{ii} = [fn{i},'{',int2str(j),'}.',fns_tmp{k}];
%                         fns{ii} = sprintf('%s{%i}.%s',fn{i},j,fns_tmp{k});
                    end
                end
            else
                ii=ii+1;
                fns{ii} = fn{i};
            end
        end
        
        % else
        %     s(is)
        %     error('allFields.m: No struct provided');
    end
end

fns = fns(1:ii);
