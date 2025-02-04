% in_matrix = arTrafoParameters(in_matrix,m,c,isdata,[trafo])
%
% Transforms derivatives with respect to the parameters when the parameter
% itself is transformed, e.g. to log10 scale
%
% Input variables:
%    in_matrix        matrix with derivatives wrt to the parameters in linear space
%    m                index of the model
%    c                index of the condition/data (depening on boolean isdata)
%    isdata           boolean, specifies if the derivatives are in a conditions 
%                     or in a data stuct. This is important since only the
%                     derivatives of the paramters in
%                     ar.model(m).condition(c).pLink / ar.model(m).data(c).pLink 
%                     are evaluated. 
%    trafo   [log10]  specifies the trafo, up to now only log10 supported
% 
% Output variables:
%    in_matrix      transformed derivatives
%
% Examples: 
%
%   syExpSimu_trafo = arTrafoParameters(ar.model(idm).data(idd).syExpSimu,idm,idd,true);
%
%       % The derivatives in ar.model(idm).data(idd).syExpSimu are wrt parameters in  
%       % linear space. The resulting syExpSimu_trafo contains the derivatives 
%       % wrt parameters in log10 space for all parameters with ar.qLog10=1.
%   

function in_matrix = arTrafoParameters(in_matrix,m,c,isdata,trafo)
global ar

dim = ndims(in_matrix);
in_size = size(in_matrix);

if(~exist('trafo','var') || isempty(trafo))
    trafo = 'log10';
end
if(~exist('isdata','var') || isempty(isdata))
    error('specify if transformation is on data or condition struct \n')
end
if(isdata)
    is_datacond = 'data';
else
    is_datacond = 'condition';
end
if(length(ar.model(m).(is_datacond)(c).pNum)>1 && in_size(end)~=length(ar.model(m).(is_datacond)(c).pNum))
    error('Check Dimensionality of your input matrix! \n')
end

%Reshape par Vector with singleton dimensions to facilitate bsxfun
pars_trafo = ar.model(m).(is_datacond)(c).pNum(ar.model(m).(is_datacond)(c).qLog10 == 1);
if(dim>2 && length(ar.model(m).(is_datacond)(c).pNum)>1)
    pars_trafo = reshape(pars_trafo,[ones(1,dim-2) size(pars_trafo)]);
end

% Assemble indexing of input matrix and
% tmp matrix that needs to be transformed   
if(length(ar.model(m).(is_datacond)(c).pNum)>1)
    which_cells = logical(false(in_size));
    sel_cells = '';
    for j=2:dim
        sel_cells = [sel_cells ':,'];
    end
    eval(['which_cells(' sel_cells 'ar.model(m).(is_datacond)(c).qLog10 == 1) = true;'])
    tmp_matrix = reshape(in_matrix(which_cells),[in_size(1:end-1) numel(pars_trafo)]);
elseif(length(ar.model(m).(is_datacond)(c).pNum)==1 && ar.model(m).(is_datacond)(c).qLog10 == 1)
    which_cells = logical(true(in_size));
    tmp_matrix = in_matrix;
elseif(length(ar.model(m).(is_datacond)(c).pNum)==1 && ar.model(m).(is_datacond)(c).qLog10 == 0)  
    return;
end

if(strcmp(trafo,'log10'))
    try
        in_matrix(which_cells) = bsxfun(@times,tmp_matrix,pars_trafo*log(10));
    catch
        tmp = tmp_matrix;
        for j=1:length(pars_trafo)
            tmp(:,:,j) = tmp(:,:,j) * pars_trafo(j)*log(10);
        end
        in_matrix(which_cells) = tmp;
    end
else
    warning('arTrafoParameters: Trying automatic transformation, cross-check or implement manual trafo! \n')
    syms f_trafo(x_trafo)
    f_trafo(x_trafo) = feval(trafo,x_trafo);
    f_diff = 1/diff(f_trafo,x_trafo);
    f_num = eval(subs(f_diff,ar.model(m).(is_datacond)(c).pNum(ar.model(m).(is_datacond)(c).qLog10 == 1)));
    if(dim>2)
        pars_trafo = reshape(f_num,[ones(1,dim-2) size(pars_trafo)]);
    end
    in_matrix(which_cells) = bsxfun(@times,tmp_matrix,pars_trafo);  
end
