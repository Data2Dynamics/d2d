% arQFit(opt,val)
% 
% This function can be used to set qFit (either for parameters or for data
% points) 
% 
%   opt     option or keyword
%           0: switch off everything
%           1: switch on everything
%           'string': The string indicates a parameter name or (if not
%           found) a name of a data set.
%           'hl': this keyword denots that data points which are
%           'highlighted' by arHighlightY
%   
%   val     value 
%           0: switch off fitting
%           1: switch on fitting
% 
% ar.qFit controls which parameters are fitted. ar.model.data.qFit controls
% which data points are used for fitting. arQfit provides some convenient
% actions to change these flags.
% 
% ar.pLabel is evaluated for the parameter names, ar.model.data.name as
% names of data sets.
% 
% Examples:
%   
% arQFit(1)
% arQFit(0)
% arQFit([],1)      % equivalent to arQFit(1)
% arQFit([],0)      % equivalent to arQFit(0)
% arQFit('data',1)  % all data sets are fitted
% arQFit('data',0)  % no data sets are fitted
% arQfit(1:10,0)
% arQfit('para_name',1)
% arQfit('data_name',1)  % change qFit of a data set, works only if parameter
%                        names does not match
% arQfit('y_name',1)     % change qFit of an observalbe, works only if parameter
%                        names and data names does not match
%
% arQfit('hl',1) set qFit=1 for the highlighted data points
% 
% See also arHighlightY


function arQFit(opt,val)
global ar

if(~exist('val','var') || isempty('val'))
    val = [];
end
if opt==0
    opt = false;
elseif opt==1
    opt = true;
end

if(isempty(opt))
    opt = 1:length(ar.qFit);
end

if(islogical(opt))
    if(length(opt)==length(ar.qFit))
        ar.qFit(find(opt)) = val;
    elseif(length(opt)==1)
        if(opt)
            ar.qFit(:)=1;
        else
            ar.qFit(:)=0;
        end
    end
elseif(isnumeric(opt))
    ar.qFit(opt) = val;
    
elseif(ischar(opt))
    done = 0;
    
    % keyword 'data'
    if(strcmp(opt,'data')==1)
        for m=1:length(ar.model)
            for d=1:length(ar.model(m).data)
                ar.model(m).data(d).qFit(:)=val;
            end
        end
        done = 1;
    end
    
    % check parameter names
    trefp = strmatch(opt,ar.pLabel,'exact');
    if(~isempty(trefp))
        ar.qFit(trefp) = val;
        disp(['qFit for parameter ',ar.pLabel{tref},' set to ',num2str(val)]);
        done = 1;
    end
    
    if(~done)  % check data sets
        for m=1:length(ar.model)
            trefd = strmatch(opt,{ar.model(m).data.name},'exact');
            for d=1:length(trefd)
                ar.model(m).data(trefd(d)).qFit(:) = val;
                disp(['qFit for data set ',ar.model(m).data(trefd(d)).name,' set to ',num2str(val)]);
                done = 1;
            end
        end
    end
    
    if(~done)  % check observables
        for m=1:length(ar.model)
            for d=1:length(ar.model(m).data)
                trefy = strmatch(opt,ar.model(m).data(d).y,'exact');
                for y=1:length(trefy)
                    ar.model(m).data(d).qFit(trefy(y)) = val;
                    disp(['qFit for observable ''',ar.model(m).data(d).y{trefy(y)}, ''' in data set ''',ar.model(m).data(d).name,''' set to ',num2str(val)]);
                    done = 1;
                end
            end
        end
    end
    
    if(~done)
        switch lower(opt)
            case {'hl','highlight','highlighted'}
                for m=1:length(ar.model)
                    for d=1:length(ar.model(m).data)
                        trefy = find(sum(ar.model(m).data(d).highlight,1)>0);
                        for y=1:length(trefy)
                            ar.model(m).data(d).qFit(trefy(y)) = val;
                            disp(['qFit for observable ''',ar.model(m).data(d).y{trefy(y)}, ''' in data set ',ar.model(m).data(d).name,''' set to ''',num2str(val)]);
                        end
                    end
                end
                
            case 'data'
                for m=1:length(ar.model)
                    for d=1:length(ar.model(m).data)
                        ar.model(m).data(d).qFit(:) = val;
                    end
                end
                
            otherwise
                error('arQFit(opt,val): opt not found.');
        end
    end
else
    error('arQFit(opt,val): opt has wrong class.');
end
