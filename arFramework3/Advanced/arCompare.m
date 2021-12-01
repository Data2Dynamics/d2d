% [same,d1,d2] = arCompare(ar1,ar2,pattern,silent)
%
%   Comparison of two ar structs/objects
%
%   ar1       an ar struct, e.g. the global variable "ar"
%   ar2       another ar struct, e.g. "ar" from an workspace
%   pattern   string or cell of strings
%             pattern is either specified by specific words
%             'Main' = {'Simu','Data','Pars','Conf','Chi2'} or 
%             evaluated by regexp.m to search for struct fields.
%             if pattern is provided, then only a subset of cells are
%             evaluated. Providing patterns is reasonable because the
%             structs often only differ by ar.fevals, ar.info, ...
%   silent    [false]
%
%   same    boolean, true if same, false if any char is different
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
% ar2.model(1).data(1).qFit = round(rand(size(ar.model(1).data(1).qFit)));
% [same,d1,d2]=arCompare(ar,ar2,{'qLog10','qFit','logfitting'})


function [same,d1,d2] = arCompare(ar1,ar2,pattern,silent)

% Load ar1,ar2 if not existent
if(nargin==0)
    filenames = fileChooserMulti('./Results', true); 
    if length(filenames)>2
       error('Error: Comparison of more than two models is not supported.') 
    end
    for j=1:length(filenames)
        fname = ['./Results/' filenames{j} '/workspace.mat'];
        if(exist(fname,'file'))
            S=load(fname);
            if j==1
                ar1 = S.ar;
            elseif j==2
                ar2 = S.ar;
            end
        else
            error('Error: No workspace found in %s',fname) 
        end    
    end
end

if(~exist('pattern','var') || isempty(pattern))
    pattern = 'all';
end
if(~exist('silent','var') || isempty(silent))
    silent = false;
end
% Original functionality
if strcmp(pattern,'all')
    [same,d1,d2] = structcmp2(ar1,ar2);
    if ~silent
        printDifference(d1,d2);
    end
% Specific comparisons for checking model differences
elseif any(strcmp(pattern,{'Main','Simu','Data','Pars','Conf','Chi2'}))
    if strcmp(pattern,'Main')
        pattern = {'Simu','Data','Pars','Conf','Chi2'};
    end
    same = 1;
    for i=1:length(pattern)
        if strcmp(pattern{i},'Simu')
            same = same * arCompareX(ar1,ar2,silent);
            same = same * arCompareZ(ar1,ar2,silent);
            same = same * arCompareV(ar1,ar2,silent);
        elseif strcmp(pattern{i},'Data')
            same = same * arCompareY(ar1,ar2,silent);
        elseif strcmp(pattern{i},'Pars')
            same = same * arComparePars(ar1,ar2);
        elseif strcmp(pattern{i},'Conf')
            same = same * arCompareConf(ar1,ar2);
        elseif strcmp(pattern{i},'Chi2')
            same = same * arCompareChi2(ar1,ar2);
        end
    end
% Specific fields in ar struct
else
    [same,d1,d2] = structcmp2(ar1,ar2,pattern);
    if ~silent
        printDifference(d1,d2);
    end    
end


function printDifference(d1,d2)
fs = allFields(d1,1);
fs = regexprep(fs,'(\w)\.(\w)','$1(1).$2');
fs = regexprep(fs,'(\w)$','$1(1)');
fs = strrep(fs,'(','{');
fs = strrep(fs,')','}');

for i=1:length(fs)
    fsname = ['.',fs{i}];  
    fsname = strrep(fsname,'{','(');
    fsname = strrep(fsname,'}',')');
    
    eval(['out1 = d1.',fs{i},';']);
    eval(['out2 = d2.',fs{i},';']);
    if isempty(out1) % do nothing
    elseif(ischar(out1))
        if(~strcmp(out1,out2))
            fprintf('%40s different: %s vs. %s \n',fsname,out1,out2);
        else
            fprintf('%40s different: %s \n',fsname,out1);
        end
    elseif(isnumeric(out1))
        fprintf('%40s different: %f vs. %f \n',fsname,out1,out2);
    end
end


% argumetn pattern can be stecified to compare only a subset of fieldnames
% pattern is compared with fieldnames via regexp.m
function [same,d1,d2] = structcmp2(s1,s2,pattern)

same = true;

% Quick check: 
% First evaluate, whether all fields are available (recursively)
fn1 = allFields(s1);
fn2 = allFields(s2);

if(exist('pattern','var') && ~isempty(pattern))
    if(~iscell(pattern))
        pattern = {pattern}; % enforce being a cell
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
    
    [sametmp,m1,m2] = valcmp(val1,val2);
    if ~sametmp 
        eval(['d1.',Fn{i},' = sprintf(''%s'',m1);']);
        eval(['d2.',Fn{i},' = sprintf(''%s'',m2);']);
    end
    same = same & sametmp; % only if both
    
%     elseif(isempty(val1) ~= isempty(val2))
%         same = false;
%         tmp = sprintf('%i x ',size(val1));
%         eval(['d1.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
%         tmp = sprintf('%i x ',size(val2));
%         eval(['d2.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
%         
%     elseif(length(size(val1)) ~= length( size(val2)))
%         same = false;
%         tmp = sprintf('%i x ',size(val1));
%         eval(['d1.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
%         tmp = sprintf('%i x ',size(val2));
%         eval(['d2.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
%         
%     elseif(sum(abs(size(val1) - size(val2) ))~=0)
%         same = false;
%         tmp = sprintf('%i x ',size(val1));
%         eval(['d1.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
%         tmp = sprintf('%i x ',size(val2));
%         eval(['d2.',Fn{i},' = ''[size = ',tmp(1:end-2),']'';']);
%         
%     elseif(isnumeric(val1) || islogical(val1))
%         if(sum(abs(val1(:)-val2(:)))>0)
%             if(length(val1)<=5)
%                 eval(['d1.',Fn{i},' = sprintf(''%f '',val1);']);
%                 eval(['d2.',Fn{i},' = sprintf(''%f '',val2);']);
%             else
%                 eval(['d1.',Fn{i},' = sprintf(''arrays different'');']);
%                 eval(['d2.',Fn{i},' = sprintf(''arrays different'');']);
%             end
%             same = false;
%         end
%         if(sum(abs(isnan(val1(:) ) - isnan( val2(:) ) ))>0)
%             same = false;
%             eval(['d1.',Fn{i},' = sprintf(''%i NaNs'',sum(isnan(val1(:))));']);
%             eval(['d2.',Fn{i},' = sprintf(''%i NaNs'',sum(isnan(val2(:))));']);
%         end
%     elseif(ischar(val1))
%         gleich = strcmp(val1,val2);
%         if(~gleich)
%             same = false;
%             val1 = val1(1:min(length(val1),100));
%             val2 = val2(1:min(length(val2),100));
%             eval(['d1.',Fn{i},' = val1;']);
%             eval(['d2.',Fn{i},' = val2;']);
%         end
%     elseif(iscell(val1))
%         if(isempty(val1) && isempty(val2))
%         else
%             try
%                 [csame,reason1,reason2] = cellcmp(val1,val2);
%             catch ERR
%                 rethrow(ERR)
%                 csame = 1;
%                 val1
%                 val2
%             end
%             if(~csame)
%                 same = false;
%                 eval(['d1.',Fn{i},' = reason1;']);
%                 eval(['d2.',Fn{i},' = reason2;']);
%             end
%         end
%     elseif(isstruct(val1))
%         error('case should not occur.')
%         %         [same2,dd1,dd2] = structcmp(s1.(fn{i}),val2);
%         %         if(~same2)
%         %             same = false;
%         %             eval(['d1.',Fn{i},' = dd1;
%         %             d2.(fn{i}) = dd2;
%         %         end
%         
%     end
end



% same = cellcmp(c1,c2)
%
%   Vergleich zwei cell Objekte. Ruft dabei valcmp auf.
function [same,reason1,reason2] = cellcmp(c1,c2)
same = true;
reason1 = '';
reason2 = '';

if isempty(c1) && isempty(c2)
    return
end
if isempty(c1) ~= isempty(c2)
    same = false;
    reason1 = sprintf('isempty=',isempty(c1));
    reason2 = sprintf('isempty=',isempty(c2));
    return
end

if iscell(c1) && iscell(c2)
    cs1 = cellfun(@class,c1,'UniformOutput',false);
    cs1 = strrep(cs1,'logical','double'); % don't distinguish these two classes
    cs2 = cellfun(@class,c2,'UniformOutput',false);
    cs2 = strrep(cs2,'logical','double'); % don't distinguish these two classes
    cscmp = strcmpCell(cs1,cs2,1);
    if sum(cscmp==0)>0
        reason1 = ['The following cells have different class: ',sprintf('%i ', find(cscmp==0))];
        reason2 = ['The following cells have different class: ',sprintf('%i ', find(cscmp==0))];
        same = false;
        return;
    end

elseif sum(abs(size(c1)-size(c2)))>0
        same = false;
        reason1 = sprintf('unequal size of cell entries');
        reason2 = sprintf('unequal size of cell entries');
        return
%     elseif size(c1,2)==1 && sum(strcmpi(cellfun(@class,c1,'UniformOutput',false)', cellfun(@class,c2,'UniformOutput',false)')~=1)>0
%         reason1 = sprintf('cellfun(@class,...) different');
%         reason2 = sprintf('cellfun(@class,...) different');
%         same = false;
else
    for i=1:numel(c1)
        [same,reason1,reason2] = valcmp(c1{i},c2{i});
        if ~same
            return;
        end
    end
end


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


function fns = allFields(s,onlyFinal)
if ~exist('onlyFinal','var') || isempty(onlyFinal)
    onlyFinal = false;
end

fns = cell(50000,1); % initialization, will be cut at the end of this function
ii=0;

for is = 1:numel(s)
    if(isstruct(s(is)))
        
        fn = fieldnames(s(is));
        for i=1:length(fn) % fieldnames, 1st level
            if(isstruct(s(is).(fn{i}))) % if again struct, do it recursively
                for j=1:length(s(is).(fn{i}))
                    fns_tmp = allFields(s(is).(fn{i})(j),onlyFinal);
                    for k=1:length(fns_tmp)
                        ii=ii+1;
                        if ~onlyFinal
                            fns{ii} = [fn{i},'(',int2str(j),').',fns_tmp{k}];
                        end
                        %                         fns{ii} = sprintf('%s(%i).%s',fn{i},j,fns_tmp{k});
                    end
                end
            elseif(iscell(s(is).(fn{i})) && ~isempty(s(is).(fn{i})))
                for j=1:length(s(is).(fn{i}))
                    if isstruct(s(is).(fn{i}){j})
                        fns_tmp = allFields(s(is).(fn{i}){j});
                        for k=1:length(fns_tmp)
                            ii=ii+1;
                            fns{ii} = [fn{i},'{',int2str(j),'}.',fns_tmp{k}];
                            %                         fns{ii} = sprintf('%s{%i}.%s',fn{i},j,fns_tmp{k});
                        end
                    else
                        ii=ii+1;
                        fns{ii} = fn{i};
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




%% Actual comparision, for cells, cellcmp is called, for structs, structcmp is called
function [same,m1,m2] = valcmp(v1,v2)
same = true;
m1 = '';
m2 = '';

cs1 = class(v1);
cs1 = strrep(cs1,'logical','double'); % don't distinguish these two classes
cs2 = class(v2);
cs2 = strrep(cs2,'logical','double'); % don't distinguish these two classes

if(strcmp(cs1,cs2)~=1)
    same = false;
    m1 = sprintf('class = %s',class(v1));
    m2 = sprintf('class = %s',class(v2));
    
elseif(isempty(v1) ~= isempty(v2))
    same = false;
    tmp = sprintf('%i x ',size(v1));
    m1 = ['size = ',tmp(1:end-2)];
    tmp = sprintf('%i x ',size(v2));
    m2 = ['size = ',tmp(1:end-2)];
    
elseif(length(size(v1)) ~= length( size(v2)))
    same = false;
    tmp = sprintf('%i x ',size(v1));
    m1 = ['size = ',tmp(1:end-2)];
    tmp = sprintf('%i x ',size(v2));
    m2 = ['size = ',tmp(1:end-2)];

elseif(sum(abs(size(v1) - size(v2) ))~=0  && sum(size(v1)==0 & size(v2)==0)==0) % no zero dims in both 
    same = false;
    tmp = sprintf('%i x ',size(v1));
    m1 = ['size = ',tmp(1:end-2)];
    tmp = sprintf('%i x ',size(v2));
    m2 = ['size = ',tmp(1:end-2)];
    
elseif isstruct(v1)
    [same,m1tmp,m2tmp] = structcmp2(v1,v2);
    if ~same
        m1 = m1tmp;
        m2 = m2tmp;
    end
    
elseif(isnumeric(v1) || islogical(v1))
    if(sum(abs(v1(:)-v2(:)))>0)
        if(length(v1)<=5)
            m1 = sprintf('%f ',v1);
            m2 = sprintf('%f ',v2);
        else
            m1 = sprintf('arrays different');
            m2 = sprintf('arrays different');
        end
        same = false;
    end
    if(sum(abs(isnan( v1(:) ) - isnan( v2(:) ) ))>0)
        same = false;
        m1 = [m1,' ',sprintf('%i NaNs',sum(isnan(v1(:))))];
        m2 = [m2,' ',sprintf('%i NaNs',sum(isnan(v2(:))))];
    end
elseif(ischar(v1))
    gleich = strcmp(v1,v2);
    if(~gleich)
        same = false;
        m1 = v1;
        m2 = v2;
    end
elseif(iscell(v1))
    if(isempty(v1) && isempty(v2))
    else
        [csame,reason1,reason2] = cellcmp(v1,v2);
        if(~csame)
            same = false;
            m1 = reason1;
            m2 = reason2;
        end
    end
elseif isa(v1,'sym')
    same = isequal(v1,v2);
    if ~same
        m1 = ['unequal symbolic expression: ',char(v1)];
        m2 = ['unequal symbolic expression: ',char(v2)];
    end    
else 
    v1
    v2
    warning('Class not yet implemented');
end
